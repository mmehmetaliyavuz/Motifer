from collections import Counter
import pandas as pd


def consensus_from_seqs(seqs, ignore_chars=set("-."), tie_char="X"):
    """Hizalanmış sekanslardan konsensüs dizi çıkarır."""
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
    MSA'dan pozisyon bazlı aa frekans matrisi (pandas DataFrame) üretir.
    Sütunlar: pos, total, A, C, D, ...
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
    Belirli eşiğin (threshold) üstünde frekansa sahip pozisyonlardan core motif çıkarır.
    """
    if freq_df.empty:
        return [], ""
    freq_df = freq_df.copy()
    # pos ve total'den sonra aa kolonları var varsayımı
    freq_df["max_freq"] = freq_df.iloc[:, 2:].max(axis=1)
    core_positions = freq_df[freq_df["max_freq"] >= threshold]["pos"].tolist()
    core_seq = "".join(consensus[p - 1] for p in core_positions)
    return core_positions, core_seq
