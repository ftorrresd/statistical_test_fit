from __future__ import annotations

from pathlib import Path


DEFAULT_INPUT_DIR = "inputs_weighted"
LEGACY_INPUT_DIR = "inputs"


def input_dir() -> Path:
    return Path(DEFAULT_INPUT_DIR)


def input_path(path: str | os.PathLike[str]) -> str:
    resolved = Path(path)
    if resolved.is_absolute():
        return str(resolved)
    if resolved.parts and resolved.parts[0] in {DEFAULT_INPUT_DIR, LEGACY_INPUT_DIR}:
        resolved = Path(*resolved.parts[1:])
    return str(input_dir() / resolved)
