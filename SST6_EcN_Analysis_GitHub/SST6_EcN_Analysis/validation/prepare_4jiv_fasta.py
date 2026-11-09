from pathlib import Path
from Bio.SeqUtils import seq1

pdb = Path("4JIV.pdb")
sequences = {chain: [] for chain in "ABCD"}

for line in pdb.read_text().splitlines():
    if line.startswith("SEQRES") and line[11] in sequences:
        sequences[line[11]].extend(line[19:].split())

chains = {chain: "".join(seq1(res) for res in residues)
          for chain, residues in sequences.items()}
assert [len(chains[c]) for c in "ABCD"] == [96, 96, 96, 93]
assert chains["A"] == chains["B"] == chains["C"]

Path("Control_4JIV.fasta").write_text(
    ">Control_4JIV_gp5x3_PAAR\n" + ":".join(chains[c] for c in "ABCD") + "\n"
)
