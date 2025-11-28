import os
import sys

# Proje kökünü PYTHONPATH'e ekle
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from motifer.pipeline import AMPAnalysisPipeline


def main():
    pipeline = AMPAnalysisPipeline(base_dir=PROJECT_ROOT)

    pipeline.download_excel()
    pipeline.load_dataframe()
    pipeline.plot_length_distribution(show=False,
                                      save_path=os.path.join(PROJECT_ROOT, "length_hist.png"))
    pipeline.export_bins_to_csv()
    pipeline.csv_to_fasta_for_bins()
    pipeline.run_mafft_for_bins()
    pipeline.compute_consensus_for_bins(threshold=0.5)

    #STREME kuruluysa aç:
    pipeline.run_streme_for_bins()

    print("[DONE] Pipeline tamamlandı.")


if __name__ == "__main__":
    main()
