from pathlib import Path
from colorir import config

PROJECT_DIR = Path(__file__).resolve().parent.parent
config.DEFAULT_PALETTES_DIR = PROJECT_DIR / "data" / "palettes"