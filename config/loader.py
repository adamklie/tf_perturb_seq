import yaml
from pathlib import Path

CONFIG_DIR = Path(__file__).parent / "colors"

def list_palettes():
    """Return available color config names (without extension)."""
    return [f.stem for f in CONFIG_DIR.glob("*.yaml")]

def load_colors(name: str = "default", subset: str = None):
    """
    Load color dictionaries for a specific dataset (or default).

    Parameters
    ----------
    name : str
        Name of the color config YAML (e.g., 'igvf_sc-islet_10X-Multiome').
    subset : str, optional
        If provided, returns only a specific subset (e.g., 'celltype_colors').

    Returns
    -------
    dict or sub-dict
        A dictionary of color mappings or a specific palette subset.
    """
    yaml_path = CONFIG_DIR / f"{name}.yaml"
    if not yaml_path.exists():
        raise FileNotFoundError(f"No color config found for '{name}' in {CONFIG_DIR}")

    with open(yaml_path) as f:
        palettes = yaml.safe_load(f)

    if subset:
        if subset not in palettes:
            raise KeyError(f"'{subset}' not found in {yaml_path.name}")
        return palettes[subset]

    return palettes
