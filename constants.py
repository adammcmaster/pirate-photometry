from pathlib import Path

OSL_DATA_PATH = Path("/STEM/data/project/osl-telescope/data/users/Pipeline/")
DATA_PATH = Path('data/')
DATA_PATH.mkdir(parents=True, exist_ok=True)

TARGET_OBSERVATIONS_PATH = DATA_PATH / 'target_observations'
TARGET_OBSERVATIONS_PATH.mkdir(exist_ok=True)