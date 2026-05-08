from __future__ import annotations
from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
from pathlib import Path
import yaml

@dataclass(frozen=True)
class ObjectiveConfig:
    magnification: float = 10.0
    NA_nominal: float = 0.25
    immersion_n: float = 1.0
    pupil_fill_ratio: float = 0.82
    focus_shift_um_850_vs_visible: float = 0.0
    defocus_um: float = 0.0

@dataclass(frozen=True)
class OCTConfig:
    center_wavelength_nm: float = 850.0
    bandwidth_nm: float = 55.0
    spectrometer_pixels: int = 2048
    source_spectrum_file: str | None = None
    pixel_to_lambda_file: str | None = None
    direct_psf_model: str = "full_spectral_rci"
    full_spectral_rci_depth_decimation: int = 1
    reference_amplitude: float = 1.0
    sample_amplitude: float = 1.0
    window: str = "hann"

@dataclass(frozen=True)
class MicroscopeConfig:
    wavelength_nm: float = 550.0
    camera_pixel_um: float = 3.45
    magnification_total: float = 10.0
    z_step_um: float = 0.25
    z_count: int = 81
    mode: str = "widefield"
    bead_diameter_um: float = 0.2
    background_adu: float = 0.0
    camera_qe_file: str | None = None
    camera_flat_field_file: str | None = None
    camera_dark_adu: float = 0.0
    camera_saturation_adu: float | None = None

@dataclass(frozen=True)
class IlluminationConfig:
    NA: float | None = None
    pupil_fill_ratio: float = 1.0
    aperture_stop_radius: float = 1.0
    field_stop_diameter_um: float | None = None

@dataclass(frozen=True)
class SampleConfig:
    scatterer_model: str = "uniform_sphere_projection"
    scatterer_diameter_um: float | None = None
    particle_refractive_index: float = 1.59
    medium_refractive_index: float = 1.33
    measurement_path: str = "epi_backscatter"
    background_amplitude: float = 0.0
    background_subtract: bool = True
    darkfield_background_fraction: float = 0.05
    scattering_lookup_file: str | None = None

@dataclass(frozen=True)
class DichroicConfig:
    angle_deg: float = 45.0
    table_file: str | None = None
    wedge_arcmin: float = 0.0
    amplitude_s: float = 1.0
    amplitude_p: float = 1.0
    phase_s_rad: float = 0.0
    phase_p_rad: float = 0.0
    mic_transmission_amplitude_s: float = 1.0
    mic_transmission_amplitude_p: float = 1.0
    mic_transmission_phase_s_rad: float = 0.0
    mic_transmission_phase_p_rad: float = 0.0
    oct_reflection_amplitude_s: float = 1.0
    oct_reflection_amplitude_p: float = 1.0
    oct_reflection_phase_s_rad: float = 0.0
    oct_reflection_phase_p_rad: float = 0.0

@dataclass(frozen=True)
class ErrorConfig:
    dichroic_tilt_deg: float = 0.0
    objective_focus_shift_um: float = 0.0
    k_linearization_rms_pixel: float = 0.0
    camera_read_noise_e: float = 2.0
    photon_gain_e_per_adu: float = 1.0
    dispersion_quadratic_rad: float = 0.0
    rolloff_per_um: float = 0.0

@dataclass(frozen=True)
class SimulationConfig:
    units: str = "um"
    random_seed: int = 20260424
    objective: ObjectiveConfig = field(default_factory=ObjectiveConfig)
    oct: OCTConfig = field(default_factory=OCTConfig)
    microscope: MicroscopeConfig = field(default_factory=MicroscopeConfig)
    illumination: IlluminationConfig = field(default_factory=IlluminationConfig)
    sample: SampleConfig = field(default_factory=SampleConfig)
    dichroic: DichroicConfig = field(default_factory=DichroicConfig)
    errors: ErrorConfig = field(default_factory=ErrorConfig)

def _merge_dataclass(cls, data: dict):
    data = data or {}
    allowed = {f.name for f in cls.__dataclass_fields__.values()}
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"Unknown keys for {cls.__name__}: {unknown}")
    return cls(**{k: v for k, v in data.items() if k in allowed})

def config_to_dict(config: SimulationConfig) -> dict:
    return asdict(config)

def config_hash(config: SimulationConfig) -> str:
    payload = json.dumps(config_to_dict(config), sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]

def _resolve_optional_path(value: str | None, config_path: Path) -> str | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return str(path.resolve())
    config_dir = config_path.parent.resolve()
    candidates = [
        config_dir / path,
        config_dir.parent / path,
        Path.cwd().resolve() / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    return str((config_dir / path).resolve())

def resolve_config_paths(config: SimulationConfig, config_path: str | Path) -> SimulationConfig:
    config_path = Path(config_path).resolve()
    return replace(
        config,
        oct=replace(
            config.oct,
            source_spectrum_file=_resolve_optional_path(config.oct.source_spectrum_file, config_path),
            pixel_to_lambda_file=_resolve_optional_path(config.oct.pixel_to_lambda_file, config_path),
        ),
        microscope=replace(
            config.microscope,
            camera_qe_file=_resolve_optional_path(config.microscope.camera_qe_file, config_path),
            camera_flat_field_file=_resolve_optional_path(config.microscope.camera_flat_field_file, config_path),
        ),
        sample=replace(
            config.sample,
            scattering_lookup_file=_resolve_optional_path(config.sample.scattering_lookup_file, config_path),
        ),
        dichroic=replace(
            config.dichroic,
            table_file=_resolve_optional_path(config.dichroic.table_file, config_path),
        ),
    )

def load_config(path: str | Path) -> SimulationConfig:
    path = Path(path).resolve()
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    config = SimulationConfig(
        units=data.get("units", "um"),
        random_seed=int(data.get("random_seed", 20260424)),
        objective=_merge_dataclass(ObjectiveConfig, data.get("objective", {})),
        oct=_merge_dataclass(OCTConfig, data.get("oct", {})),
        microscope=_merge_dataclass(MicroscopeConfig, data.get("microscope", {})),
        illumination=_merge_dataclass(IlluminationConfig, data.get("illumination", {})),
        sample=_merge_dataclass(SampleConfig, data.get("sample", {})),
        dichroic=_merge_dataclass(DichroicConfig, data.get("dichroic", {})),
        errors=_merge_dataclass(ErrorConfig, data.get("errors", {})),
    )
    return resolve_config_paths(config, path)

def dump_resolved_config(path: str | Path, config: SimulationConfig) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(config_to_dict(config), sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
