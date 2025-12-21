import os
import sys

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


    results = pipeline.compute_consensus_for_bins(threshold=0.5)

    pipeline.export_consensus(
        results,
        os.path.join(PROJECT_ROOT, "consensus_and_core.fasta")
    )


    out_txt = os.path.join(PROJECT_ROOT, "core_motifs_threshold_0.5.txt")
    with open(out_txt, "w") as f:
        for (low, high), d in results.items():
            f.write(f"{low}-{high} aa\n")
            f.write(f"core_positions: {d['core_positions']}\n")
            f.write(f"core_sequence : {d['core_sequence']}\n")
            f.write(f"consensus     : {d['consensus']}\n")
            f.write("\n")

    print(f"[OK] Core motifs are wrote: {out_txt}")

 

    print("[DONE] Pipeline is done.")



if __name__ == "__main__":
    main()

