from collections import Counter
import pandas as pd


def consensus_from_seqs(seqs, ignore_chars=set("-."), tie_char="X"):
    """It extracts a consensus sequence from aligned sequences."""
    if not seqs:
        return ""
    L = max(len(s) for s in seqs)
    consensus = []
    for i in range(L):
        col = []
        for s in seqs:
            if i < len(s):
                aa = s[i]
                if aa not in ignore_chars:
                    col.append(aa)
        if not col:
            consensus.append("-")
            continue
        counts = Counter(col)
        mc = counts.most_common()
        if len(mc) > 1 and mc[0][1] == mc[1][1]:
            consensus.append(tie_char)
        else:
            consensus.append(mc[0][0])
    return "".join(consensus)


def column_frequency_matrix(
    seqs,
    alphabet="ACDEFGHIKLMNPQRSTVWY",
    ignore_chars="-.",
):
    """
   Generates a position-based frequency matrix (pandas DataFrame) from MSA. Columns: pos, total, A, C, D, ...
    """
    if not seqs:
        return pd.DataFrame()
    L = max(len(s) for s in seqs)
    letters = list(alphabet)
    rows = []

    for pos in range(L):
        col = []
        for s in seqs:
            if pos < len(s):
                aa = s[pos]
                if aa not in ignore_chars:
                    col.append(aa)
        total = len(col)
        counts = Counter(col)
        row = {"pos": pos + 1, "total": total}
        for aa in letters:
            row[aa] = counts[aa] / total if total > 0 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def extract_motif_core(consensus, freq_df, threshold=0.50):
    """
It extracts core motifs from positions with frequencies above a certain threshold.
    """
    if freq_df.empty:
        return [], ""
    freq_df = freq_df.copy()
    # Assumption that there are columns aa after pos and total
    freq_df["max_freq"] = freq_df.iloc[:, 2:].max(axis=1)
    core_positions = freq_df[freq_df["max_freq"] >= threshold]["pos"].tolist()
    core_seq = "".join(consensus[p - 1] for p in core_positions)
    return core_positions, core_seq
