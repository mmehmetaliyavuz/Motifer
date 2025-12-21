from collections import Counter

def normalize_seq(seq: str) -> str:
    """Normalizing simple amino acids; U->C, O->K, etc."""
    seq = str(seq).upper().replace(" ", "")
    seq = seq.replace("U", "C").replace("O", "K")
    return seq


def read_fasta_alignment(path):
    """The sequence list is returned from the FASTA alignment file."""
    seqs = []
    current = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
            else:
                current.append(line)
        if current:
            seqs.append("".join(current))
    return seqs
