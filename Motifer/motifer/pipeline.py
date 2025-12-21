import os
import math
import subprocess

import pandas as pd
import requests
import matplotlib.pyplot as plt

from .config import BASE_DIR, DRAMP_EXCEL_URL, DRAMP_EXCEL_FILENAME, LENGTH_BINS
from .io_utils import normalize_seq, read_fasta_alignment
from .motif_utils import consensus_from_seqs, column_frequency_matrix, extract_motif_core


class AMPAnalysisPipeline:
    def __init__(self, base_dir: str | None = None, length_bins=None):
        self.base_dir = base_dir or BASE_DIR
        self.length_bins = length_bins or LENGTH_BINS
        self.df: pd.DataFrame | None = None

    # ---------- 1. Data ----------

    def download_excel(self):
        path = os.path.join(self.base_dir, DRAMP_EXCEL_FILENAME)
        if os.path.exists(path):
            print(f"[INFO] Excel is already exist: {path}")
            return path

        print(f"[INFO] DRAMP Excel is downloading: {DRAMP_EXCEL_URL}")
        resp = requests.get(DRAMP_EXCEL_URL)
        if resp.status_code == 200:
            with open(path, "wb") as f:
                f.write(resp.content)
            print(f"[OK] Recorded: {path}")
        else:
            raise RuntimeError(f"Excel could not dowloaded (status={resp.status_code})")
        return path

    def load_dataframe(self):
        excel_path = os.path.join(self.base_dir, DRAMP_EXCEL_FILENAME)
        print(f"[INFO] Excel is reading: {excel_path}")
        df = pd.read_excel(excel_path)

        if "Sequence" not in df.columns:
            raise KeyError("The DataFrame does not have a 'Sequence' column.")

        df = df.dropna(subset=["Sequence"]).copy()
        df["Sequence"] = df["Sequence"].astype(str)
        df["length"] = df["Sequence"].str.len()
        self.df = df
        print(f"[OK] Total peptide: {len(df)}")
        return df

    # ---------- 2. Lenght Range  ----------

    def plot_length_distribution(self, bin_size=3, show=True, save_path=None):
        if self.df is None:
            raise ValueError("There is no data frame. Call load_dataframe() first..")

        max_len = self.df["length"].max()
        max_bin = int(math.ceil(max_len / bin_size) * bin_size)
        bins = list(range(0, max_bin + bin_size, bin_size))
        labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1)]

        self.df["length_bin"] = pd.cut(
            self.df["length"],
            bins=bins,
            labels=labels,
            include_lowest=True,
            right=True,
        )
        counts = self.df["length_bin"].value_counts().sort_index()

        plt.figure(figsize=(10, 5))
        plt.bar(counts.index.astype(str), counts.values)
        plt.xlabel("Lenght Range (aa)")
        plt.ylabel("AMP Number")
        plt.title("AMP Length Distribution")
        plt.xticks(rotation=45)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"[OK] Histogram recorded: {save_path}")
        if show:
            plt.show()

    # ---------- 3. Binlere göre CSV ----------

    def export_bins_to_csv(self):
        if self.df is None:
            raise ValueError("There is no data frame. Call load_dataframe() first..")

        for low, high in self.length_bins:
            subset = self.df[
                (self.df["length"] >= low) & (self.df["length"] <= high)
            ].copy()
            fname = os.path.join(self.base_dir, f"{low}_{high}aa_peptides.csv")
            subset.to_csv(fname, index=False)
            print(f"[OK] {fname} are created → {len(subset)} peptide")

    # ---------- 4. CSV -> FASTA ----------

    def csv_to_fasta(self, csv_path, fasta_path, seq_col="Sequence", id_col=None):
        df = pd.read_csv(csv_path)
        df = df.dropna(subset=[seq_col]).copy()

        with open(fasta_path, "w") as f:
            for i, row in df.iterrows():
                raw_seq = row[seq_col]
                seq = normalize_seq(raw_seq)

                if id_col and id_col in df.columns and pd.notna(row[id_col]):
                    header = str(row[id_col])
                else:
                    header = f"seq{i+1}_len{len(seq)}"

                f.write(f">{header}\n{seq}\n")
        print(f"[OK] FASTA oluşturuldu: {fasta_path}")

    def csv_to_fasta_for_bins(self):
        for low, high in self.length_bins:
            csv_file = os.path.join(self.base_dir, f"{low}_{high}aa_peptides.csv")
            fasta_file = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_normalized.fasta"
            )
            self.csv_to_fasta(csv_file, fasta_file)

    # ---------- 5. MAFFT ----------

    def run_mafft(self, input_fa, output_fa, anysymbol=False):
        cmd = ["mafft", "--auto"]
        if anysymbol:
            cmd.insert(1, "--anysymbol")
        cmd.append(input_fa)

        print(f"[CMD] {' '.join(cmd)} > {output_fa}")
        with open(output_fa, "w") as fout:
            subprocess.run(cmd, stdout=fout, check=True)

    def run_mafft_for_bins(self):
        for low, high in self.length_bins:
            in_fa = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_normalized.fasta"
            )
            out_fa = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_mafft.fasta"
            )
            print(f"[INFO] MAFFT: {in_fa} → {out_fa}")
            self.run_mafft(in_fa, out_fa, anysymbol=False)

    # ---------- 6. Consensüs & core motif ----------

    def compute_consensus_for_bins(self, threshold=0.5):
        results = {}
        for low, high in self.length_bins:
            aln_path = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_mafft.fasta"
            )
            print(f"\n=== {low}-{high} aa ===")
            seqs = read_fasta_alignment(aln_path)
            print(f"  {len(seqs)} sekans.")

            consensus = consensus_from_seqs(seqs)
            freq_df = column_frequency_matrix(seqs)
            core_pos, core_seq = extract_motif_core(consensus, freq_df, threshold)

            print("  Core position (fisrt 20):", core_pos[:20])
            print("  Core motif (fisrt 50 aa):", core_seq[:50])

            results[(low, high)] = {
                "consensus": consensus,
                "freq_df": freq_df,
                "core_positions": core_pos,
                "core_sequence": core_seq,
            }
        return results

    def export_consensus(self, results, out_fasta):
        with open(out_fasta, "w") as f:
            for (low, high), d in results.items():
                f.write(f">consensus_{low}_{high}\n{d['consensus']}\n")
                f.write(f">core_{low}_{high}\n{d['core_sequence']}\n")

    # ---------- 7. STREME  ----------

    def run_streme_for_bins(
        self,
        minw=3,
        maxw=8,
        evalue=1e-3,
        nmotifs=5,
        seed=42,
    ):
        """
       This assumes that STREME is installed in the system PATH. This method will only work if the `streme` command exists.

        """
        import shutil

        if shutil.which("streme") is None:
            raise RuntimeError(
                "'The 'streme' command was not found. Make sure MEME Suite is installed and in the PATH."
            )

        for low, high in self.length_bins:
            fasta_in = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_normalized.fasta"
            )
            fasta_clean = os.path.join(
                self.base_dir, f"{low}_{high}aa_peptides_streme.fasta"
            )
            out_dir = os.path.join(
                self.base_dir, f"streme_{low}_{high}_out"
            )

            valid_aas = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
            total = 0
            cleaned = 0
            with open(fasta_in, "r") as fin, open(fasta_clean, "w") as fout:
                for line in fin:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        fout.write(line + "\n")
                        continue
                    if not line.strip():
                        continue
                    total += 1
                    seq = line.strip().upper()
                    new_seq = "".join(
                        ch if ch in valid_aas else "X" for ch in seq
                    )
                    if new_seq != seq:
                        cleaned += 1
                    fout.write(new_seq + "\n")
            print(
                f"[INFO] {low}-{high}aa cleaned FASTA: {fasta_clean} "
                f"(total {total}, cleaned  {cleaned})"
            )

            cmd = [
                "streme",
                "-p",
                fasta_clean,
                "-oc",
                out_dir,
                "-minw",
                str(minw),
                "-maxw",
                str(maxw),
                "--protein",
                "--evalue",
                "--thresh",
                str(evalue),
                "--nmotifs",
                str(nmotifs),
                "--seed",
                str(seed),
            ]
            print(f"[CMD] {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            print(f"[OK] STREME results: {out_dir}")
