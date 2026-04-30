from __future__ import annotations
from dataclasses import replace
from .config_schema import SimulationConfig, ErrorConfig

def with_error(config: SimulationConfig, **kwargs) -> SimulationConfig:
    err = replace(config.errors, **kwargs)
    return replace(config, errors=err)

def make_one_factor_sweep(config: SimulationConfig, name: str, values):
    for v in values:
        yield with_error(config, **{name: v})
