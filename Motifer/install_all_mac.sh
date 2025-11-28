#!/usr/bin/env bash
set -e

echo "========================================"
echo "  MOTIFER FULL INSTALLER (macOS)"
echo "  Python venv + pip deps + MAFFT + MEME"
echo "========================================"

# ---- CONFIG ----
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$PROJECT_ROOT"

# ---- CHECK HOMEBREW ----
if ! command -v brew >/dev/null 2>&1; then
  echo "[INFO] Homebrew bulunamadı. Kuruluyor..."
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
else
  echo "[OK] Homebrew bulundu."
fi

echo "[INFO] Homebrew güncelleniyor..."
brew update

# ---- INSTALL MAFFT ----
if ! command -v mafft >/dev/null 2>&1; then
  echo "[INFO] MAFFT kuruluyor..."
  brew install mafft
else
  echo "[OK] MAFFT zaten kurulu."
fi

# ---- INSTALL MEME SUITE (includes STREME) ----
# ---- MEME SUITE / STREME (MANUEL) ----
if ! command -v streme >/dev/null 2>&1; then
  echo "[WARN] 'streme' komutu bulunamadı."
  echo "       Homebrew'de 'meme' formülü yok, o yüzden otomatik kurulamıyor."
  echo "       MEME Suite'i resmi siteden manuel kurman gerekiyor:"
  echo "       https://meme-suite.org/meme/"
else
  echo "[OK] MEME Suite + STREME zaten kurulu."
fi


# ---- PYTHON ENV ----
echo "[INFO] Python virtual environment oluşturuluyor (.venv)..."
python3 -m venv .venv

echo "[INFO] venv aktif ediliyor..."
source .venv/bin/activate

echo "[INFO] pip yükseltiliyor..."
pip install --upgrade pip

echo "[INFO] Python bağımlılıkları yükleniyor..."
if [ -f "requirements.txt" ]; then
  pip install -r requirements.txt
else
  echo "[WARN] requirements.txt bulunamadı, atlanıyor."
fi

echo "========================================"
echo "    🎉 KURULUM TAMAMLANDI 🎉"
echo "========================================"
echo ""
echo "Kullanım:"
echo "  source .venv/bin/activate"
echo "  python scripts/run_full_pipeline.py"
echo ""
echo "MAFFT path:"
which mafft
echo "STREME path:"
which streme


