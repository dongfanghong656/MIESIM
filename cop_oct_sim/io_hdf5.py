from __future__ import annotations
from pathlib import Path
import h5py
import json
import numpy as np

def _write_item(group, key: str, value):
    if isinstance(value, dict):
        child = group.create_group(key)
        for sub_key, sub_value in value.items():
            _write_item(child, str(sub_key), sub_value)
        return
    if isinstance(value, str):
        group.attrs[key] = value
        return
    if value is None:
        group.attrs[key] = "None"
        return
    arr = np.asarray(value)
    kwargs = {}
    if arr.ndim > 0 and arr.size > 64:
        kwargs = {"compression": "gzip", "compression_opts": 4, "shuffle": True}
    group.create_dataset(key, data=arr, **kwargs)

def write_dict_h5(path, data: dict, metadata: dict | None = None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for k, v in data.items():
            _write_item(f, str(k), v)
        if metadata:
            f.attrs["metadata_json"] = json.dumps(metadata, ensure_ascii=False)
