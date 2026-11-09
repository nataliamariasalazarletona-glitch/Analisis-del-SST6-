from Bio.PDB import PDBParser, Superimposer
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
from itertools import permutations

parser = PDBParser(QUIET=True)

exp = parser.get_structure(
    "experimental", "4JIV_experimental.pdb"
)[0]

pred = parser.get_structure(
    "predicted", "4JIV_predicted_rank001.pdb"
)[0]

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

def get_residues(chain):
    return [
        r for r in chain
        if r.id[0] == " " and "CA" in r
    ]

def get_sequence(residues):
    return "".join(seq1(r.resname) for r in residues)

def matched_atoms(exp_chain, pred_chain):
    er = get_residues(exp_chain)
    pr = get_residues(pred_chain)

    alignment = aligner.align(
        get_sequence(er), get_sequence(pr)
    )[0]

    exp_atoms = []
    pred_atoms = []

    for (e0, e1), (p0, p1) in zip(*alignment.aligned):
        for i, j in zip(range(e0, e1), range(p0, p1)):
            if er[i].resname == pr[j].resname:
                exp_atoms.append(er[i]["CA"])
                pred_atoms.append(pr[j]["CA"])

    return exp_atoms, pred_atoms

results = []

for perm in permutations("ABC"):
    exp_atoms = []
    pred_atoms = []

    for ec, pc in zip("ABC", perm):
        ea, pa = matched_atoms(exp[ec], pred[pc])
        exp_atoms.extend(ea)
        pred_atoms.extend(pa)

    ea, pa = matched_atoms(exp["D"], pred["D"])
    exp_atoms.extend(ea)
    pred_atoms.extend(pa)

    sup = Superimposer()
    sup.set_atoms(exp_atoms, pred_atoms)

    results.append((perm, sup.rms, len(exp_atoms)))

results.sort(key=lambda x: x[1])

print("====================================")
print("VALIDACION ESTRUCTURAL PDB 4JIV")
print("====================================")

for perm, rmsd, n in results:
    print(
        f"Permutacion gp5: {perm} | "
        f"RMSD: {rmsd:.3f} A | "
        f"Atomos CA: {n}"
    )

best = results[0]

print("\nMEJOR RESULTADO")
print("Permutacion:", best[0])
print("RMSD:", round(best[1], 3), "A")
print("Atomos comparados:", best[2])
