# Common-Path OCT Simulation

二向色镜显微镜-共光路 SD-OCT 的第一阶段仿真项目。

核心边界：

1. `microscope_forward` 独立生成显微镜理论/实测 PSF；
2. `oct_forward` 独立生成 OCT raw interferogram 和 direct OCT PSF；
3. 显微镜到 OCT 的转换公式只放在 `theory_conversion`，作为被测对象；
4. `validation` 和 `scripts/run_validation_suite.py` 生成测试数据包和指标表。

快速运行：

```powershell
python -m pytest
python scripts/run_validation_suite.py --output-root outputs
python scripts/run_validation_suite.py --config configs/config_full_spectral_rci_smoke.yaml --output-root outputs --N 12 --k-samples 24 --pad-factor 1
```

验证输出会写入 `outputs/validation_YYYYMMDD_HHMMSS`。

当前验证器会生成：

- direct OCT raw / PSF HDF5；
- low-resolution full spectral RCI direct PSF HDF5，用于审查 `h_RCI(x,y,z;k)` 到 k-space FFT 的直接链路；
- OCT reconstruction HDF5，包含 `complex`、`magnitude`、`depth_um`、`k_linear`；
- 显微镜 z-stack HDF5；
- 路径 E measured predictive PSF 与 ideal upper-bound PSF；
- `figures/*.png`，包含显微 XY/XZ/YZ、OCT A-scan、Path E residual 和 sweep trend；
- `validation_summary.json`；
- 覆盖 pupil fill、finite bead、focus shift、dispersion、k-linearization、roll-off 的 `sweep_summary.csv`。
- `negative_controls.csv`，故意注入 chromatic mismatch、wrong dichroic path、bad k calibration 并确认会被检测出来；
- `convergence_summary.csv`，覆盖 N、pad_factor、k_samples、z_step、window、bead size、camera pixel 的收敛/敏感性扫描；
- `direct_model_comparison.csv`，低分辨率比较 `hybrid_rci` 与 `full_spectral_rci` direct PSF，作为 full spectral RCI 路径的 diagnostic evidence；
- `direct_model_axial_profiles.csv`，逐 depth 输出 `full_spectral_rci` 与 `hybrid_rci` 的 max / on-axis / integrated axial profiles 和残差；
- machine-readable `verdict`，包含 `pilot_pass`、`review_required` 和 `unsupported_claims`。

当前 gate 语义：

- `checks.all_pass` 只表示所有已编码 gate 通过；
- `verdict.pilot_pass` 只有在 `checks.all_pass=True` 且 `blocker_count=0` 时才为 true；
- `verdict.review_required=True` 表示仍有 blocker 或 unsupported claim，不得解释为物理真值已验证；
- `direct_model_comparison.csv` 中 hybrid/full spectral axial 或 lateral FWHM 相对误差超过阈值时，`direct_model_comparison.pass=False`。

数据产品中的关键审查字段：

- `microscope_zstack.h5`：`ideal`、`raw`、`background`、`corrected`、`bead_convolved`、`camera_integrated`、`measured`、`x_um`、`y_um`、`z_um`；
- `microscope_zstack.h5` 还会记录 `scatterer_model`、`measurement_path`、`scattering_raw`、`scattering_background`、`scattering_corrected`、`finite_target_fwhm_ratio`、`camera_qe`、`camera_flat_field`、`camera_dark_adu`、`camera_saturation_adu`、`illumination_power_fraction`；
- `microscope_forward/illumination_pupil.h5`：`illumination_pupil`、`illumination_psf`、`illumination_envelope`、`illumination_na`、`illumination_pupil_fill_ratio`、`illumination_aperture_stop_radius`；
- `oct_raw_interferogram.h5`：`pixel_index`、`wavelength_nm`、`k_true`、`k`、`k_linear`、`source`、`reference_field`、`sample_field`、`background`、`interference`、`raw`、`h_rci_sample`、`sample_scattering_amplitude`；
- `oct_psf_direct.h5`：`psf`、`magnitude`、`rci_through_focus`、`separability_nrmse`、`full_spectral_rci`、`depth_um`、`x_um`、`y_um`、`k_true`、`k_measured`、`k_linear`。
- `oct_psf_direct_lowres_full_spectral_rci.h5`：低分辨率 `spectral_rci_direct_psf` 审查产物；同时记录未加 source/window/dispersion/fixed scatterer-depth phase 前的 `spectral_rci_sampled_h_rci[depth,k,y,x]`、`spectral_rci_k_true`、`spectral_rci_k_measured`、`spectral_rci_k_linear`、`spectral_rci_wavelength_nm`、`spectral_rci_dx_um_per_k`、`spectral_rci_source`、`spectral_rci_window`、`spectral_rci_dispersion_phase_rad`、`spectral_rci_depth_phase_origin_um`、`spectral_rci_axial_profile_max`、`spectral_rci_axial_profile_on_axis`、`spectral_rci_axial_profile_integrated`；当启用 depth decimation 时还会记录 `spectral_rci_computed_depth_indices`、`spectral_rci_computed_depth_um`、`full_spectral_rci_depth_decimation` 和 `full_spectral_rci_interpolated`。模型约定、phase sign 和 normalization 作为 HDF5 attributes 写入。

显微镜相机模型会先在物方坐标中进行有限球体卷积，再按 camera pixel 的物方尺寸做 block integration；当相机像素大于仿真采样间隔时，`measured` 会落在降采样后的 `camera_x_um` / `camera_y_um` 网格上。随后会应用 camera QE、flat field、background/dark frame、shot/read noise 和可选 saturation。

照明模型由 `illumination` 配置段控制。`NA` 和 `aperture_stop_radius` 限制照明 pupil support，`pupil_fill_ratio` 控制照明光瞳填充，`field_stop_diameter_um` 生成样品平面视场包络；当照明 NA、孔径或填充显式受限时，显微 forward 会把 illumination PSF 作用到 epi scatterer response 的形状上，而不只是缩放强度。validation 会单独输出 `illumination_pupil.h5` 供审查。

散射球样品模型由 `sample` 配置段控制：

- `scatterer_model` 支持 `delta`、`gaussian`、`uniform_sphere_projection`、`rayleigh`、`rayleigh_mie_lookup`；
- `measurement_path` 支持 `epi_backscatter` 与 `oblique_darkfield`；
- `rayleigh_mie_lookup` 必须提供 `scattering_lookup_file`，CSV 至少包含 `wavelength_nm,amplitude`，可选 `phase_rad`；
- `finite_target_fwhm_ratio` 会报告球直径相对显微中心平面 FWHM 的比例，用于检查散射球是否仍可近似点目标。

direct OCT PSF 支持两种模型，由 `oct.direct_psf_model` 控制：

- `hybrid_rci`：默认快速模型，生成 `rci_through_focus[z,y,x]`，再与谱域重建得到的轴向 gate 组合；适合常规 validation 和参数扫描。
- `full_spectral_rci`：低分辨率全光谱模型，把 `h_RCI(x,y,z;k)` 放在 source spectrum、window、dispersion phase、depth phase、k-linearization 和 k-space FFT 前生成 direct truth 审查链。该路径更接近审查要求，但目前仍是 scalar low-NA、低分辨率实现，尚不是生产级全物理真值。

`separability_nrmse` 用于量化可分离近似的残差。

`oct.full_spectral_rci_depth_decimation` 可用于扩展 full spectral RCI 计算：`1` 表示每个 depth plane 都直接计算；大于 `1` 时只直接计算抽样深度并插值回完整 depth grid，同时保留首点、中点、末点作为 anchor。

defocus 约定：

- `objective.defocus_um` 是 system wavefront Zernike defocus OPD coefficient，不是样品物理 z 位移；
- 样品/焦点物理 z 位移通过 low-NA physical defocus phase 进入 propagation；
- `errors.objective_focus_shift_um` 和 `objective.focus_shift_um_850_vs_visible` 表示 physical focus shift，不会再作为 Zernike OPD coefficient 使用。

配置文件中的 CSV/YAML 相对路径会在 `load_config(path)` 中解析为绝对路径；解析顺序优先尝试 config 文件所在目录，其次尝试项目根目录，因此 `configs/config_visible_660.yaml` 可从项目根目录或子目录运行。

若 `SimulationConfig.oct.pixel_to_lambda_file` 或 `SimulationConfig.oct.source_spectrum_file` 指向 CSV 文件，OCT direct forward 会优先使用实测像素-波长标定和光源谱；否则使用配置中的高斯谱低维模型。

若 `SimulationConfig.microscope.camera_qe_file` 指向 CSV 文件，显微镜 forward 会按 `wavelength_nm,qe` 或 `wavelength_nm,response` 插值相机量子效率，并把结果写入 `camera_qe`。

Path E 现在分为两类输出：

- `converters/path_e_measured_predictive_psf.h5`：空间项来自显微 measured/corrected stack，轴向项由 source spectrum、k grid、dispersion 和 reconstruction contract 独立预测；
- `converters/path_e_ideal_upper_bound_psf.h5`：空间项和轴向项使用理想/直接信息，只作为上限/自洽对照；
- legacy aliases 仍输出 `path_e_predictive_psf.h5` 与 `path_e_oracle_upper_bound_psf.h5`，但审查时应优先看新名称。

Path E pass gate 不再只看 3D NRMSE，还要求 axial 与 lateral FWHM 一致；因此低 NRMSE 但轴向 FWHM 错误的情况会被判为 review。

二向色镜模型已区分 `mic_transmission` 和 `oct_reflection`。若 `dichroic.table_file` 指向 CSV 文件，可用 `wavelength_nm, AOI, Rs, Rp, Ts, Tp, phase_s, phase_p` 表驱动显微透射路径和 OCT 反射路径。
