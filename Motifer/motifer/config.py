import os

# Varsayılan temel dizin (isteğe göre değiştir)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DRAMP_EXCEL_URL = (
    "http://dramp.cpu-bioinfor.org/downloads/download.php"
    "?filename=download_data/DRAMP3.0_new/general_amps.xlsx"
)

DRAMP_EXCEL_FILENAME = "general_amps.xlsx"

# Uzunluk binleri
LENGTH_BINS = [
    (9, 15),
    (12, 18),
    (15, 25),
    (20, 35),
    (25, 40),
]
