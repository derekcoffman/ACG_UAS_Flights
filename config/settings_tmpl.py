import os

BASE_DATA_PATH = "/base/path/to/data"

ENVDSYS_DATA_PATH = os.path.join(BASE_DATA_PATH, "envDataSystem/from_cloudbase")

ENVDSYS_CLOUDYSKY_DATA_PATH = os.path.join(ENVDSYS_DATA_PATH, "UIServer/cloudy1.acg.pmel.noaa.gov/cloudysky")
ENVDSYS_GROUDSTATION_DATA_PATH = ENVDSYS_DATA_PATH

data_paths = {
    "base_data_path": BASE_DATA_PATH,
    "envdsys_data_path": ENVDSYS_DATA_PATH,
    "envdsys_cloudysky_data_path": ENVDSYS_CLOUDYSKY_DATA_PATH,
    "envdsys_groundstation_data_path": ENVDSYS_GROUDSTATION_DATA_PATH
}