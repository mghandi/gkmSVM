#!/usr/bin/env python3
"""Generate the golden-test fixture FASTAs deterministically (seeded PRNG, no downloads).

The generated files are committed; re-running this script must reproduce them byte-for-byte
(`python3 tests/golden/make_fixtures.py --check`). Fixtures are deliberately small so that the frozen
expected outputs stay small, and each file exercises one edge of the FASTA reader or the algorithm:

  dna_small_{pos,neg}.fa   12+12 x 100 bp DNA, a 7-mer motif planted in every positive   (main grid)
  dna_medium_{pos,neg}.fa  40+40 x 300 bp DNA                                             (threading, alg=1)
  oracle_{pos,neg}.fa      6+6 x 60 bp DNA                                                (exact-arithmetic oracle)
  dupnames_pos.fa          12 records, records 0 and 1 share the name "posX"              (silent merge, issue 1)
  longnames_pos.fa         12 records whose names are 120 characters long                 (name buffer overflow)
  emptyrec_pos.fa          12 records; #3 has an empty sequence, #7 is all 'N'            (records skipped)
  lowercase_pos.fa         dna_small_pos.fa in lower case                                 (case folding)
  withN_pos.fa             dna_small_pos.fa with 'N' runs spliced in                      (non-alphabet chars dropped)
  crlf_pos.fa              dna_small_pos.fa with CRLF line endings
  multiline_pos.fa         dna_small_pos.fa wrapped at 40 columns
  b2_{pos,neg}.fa          10+10 x 80 over the 2-letter alphabet "AB" (alphabet_b2.txt)
  b5_{pos,neg}.fa          6+6 x 60 over "ACGTN" as a 5-letter alphabet (alphabet_b5.txt)  (crashes today; Phase 5)
  test_seqs.fa             10 x 100 bp DNA test sequences for gkmsvm_classify
  sv_svseq.fa / sv_svalpha.out  a hand-made support-vector model: 5 pos + 5 neg from dna_small with fixed alphas
  sv_longnames_* / test_longnames.fa  the same model and 10 test sequences with 120-character names   (Phase 1)
  sv_b2_* / test_b2.fa     a fixed model and 8 test sequences over the "AB" alphabet                  (Phase 1)
  names_edge_pos.fa        empty / tab / non-ASCII / 200-char / cross-file-duplicate names              (Phase 3)
  sv_dupalpha_svalpha.out, sv_missing_svalpha.out  invalid SV models (duplicate / missing name)         (Phase 3)
  protein_{pos,neg}.fa + alphabet_protein.txt   6+6 x 60 over the 20 amino acids                          (Phase 5)
  rna_{pos,neg}.fa + alphabet_rna_pairs.txt     dna_small with T->U, alphabet file declaring complements  (Phase 5)
  alphabet_b5_pairs_bad.txt                     complement declared for some symbols only (rejected)      (Phase 5)
"""
import argparse, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
MOTIF = "TGACTCA"


def seqs(rng, n, length, alpha="ACGT", motif=None):
    out = []
    for _ in range(n):
        s = [rng.choice(alpha) for _ in range(length)]
        if motif:
            p = rng.randrange(0, length - len(motif))
            s[p:p + len(motif)] = list(motif)
        out.append("".join(s))
    return out


def fasta(records, eol="\n", width=None):
    parts = []
    for name, s in records:
        parts.append(">" + name + eol)
        if width:
            for i in range(0, len(s), width):
                parts.append(s[i:i + width] + eol)
        else:
            parts.append(s + eol)
    return "".join(parts)


def build():
    rng = random.Random(20260829)
    files = {}
    pos = seqs(rng, 12, 100, motif=MOTIF)
    neg = seqs(rng, 12, 100)
    small_pos = [(f"pos{i}", s) for i, s in enumerate(pos)]
    small_neg = [(f"neg{i}", s) for i, s in enumerate(neg)]
    files["dna_small_pos.fa"] = fasta(small_pos)
    files["dna_small_neg.fa"] = fasta(small_neg)
    files["dna_medium_pos.fa"] = fasta([(f"mp{i}", s) for i, s in enumerate(seqs(rng, 40, 300, motif=MOTIF))])
    files["dna_medium_neg.fa"] = fasta([(f"mn{i}", s) for i, s in enumerate(seqs(rng, 40, 300))])
    files["oracle_pos.fa"] = fasta([(f"op{i}", s) for i, s in enumerate(seqs(rng, 6, 60, motif=MOTIF))])
    files["oracle_neg.fa"] = fasta([(f"on{i}", s) for i, s in enumerate(seqs(rng, 6, 60))])
    files["dupnames_pos.fa"] = fasta([("posX" if i < 2 else n, s) for i, (n, s) in enumerate(small_pos)])
    files["longnames_pos.fa"] = fasta([("L" * 116 + f"_{i:03d}", s) for i, (n, s) in enumerate(small_pos)])
    emp = list(small_pos)
    emp[3] = ("pos3_empty", "")
    emp[7] = ("pos7_allN", "N" * 100)
    files["emptyrec_pos.fa"] = fasta(emp)
    files["lowercase_pos.fa"] = fasta([(n, s.lower()) for n, s in small_pos])
    withn = []
    for i, (n, s) in enumerate(small_pos):
        p = 10 + 5 * i
        withn.append((n, s[:p] + "NNN" + s[p:p + 30] + "N" + s[p + 30:]))
    files["withN_pos.fa"] = fasta(withn)
    files["crlf_pos.fa"] = fasta(small_pos, eol="\r\n")
    files["multiline_pos.fa"] = fasta(small_pos, width=40)
    files["b2_pos.fa"] = fasta([(f"b2p{i}", s) for i, s in enumerate(seqs(rng, 10, 80, alpha="AB", motif="ABBABAA"))])
    files["b2_neg.fa"] = fasta([(f"b2n{i}", s) for i, s in enumerate(seqs(rng, 10, 80, alpha="AB"))])
    files["alphabet_b2.txt"] = "A\nB\n"
    files["b5_pos.fa"] = fasta([(f"b5p{i}", s) for i, s in enumerate(seqs(rng, 6, 60, alpha="ACGTN", motif=MOTIF))])
    files["b5_neg.fa"] = fasta([(f"b5n{i}", s) for i, s in enumerate(seqs(rng, 6, 60, alpha="ACGTN"))])
    files["alphabet_b5.txt"] = "A\nC\nG\nT\nN\n"
    files["test_seqs.fa"] = fasta([(f"t{i}", s) for i, s in enumerate(seqs(rng, 10, 100, motif=MOTIF if rng.random() < 0.5 else None))])
    # a fixed (not trained) support-vector model in the {prefix}_svalpha.out / {prefix}_svseq.fa format
    sv = [small_pos[i] for i in (0, 2, 5, 7, 11)] + [small_neg[i] for i in (1, 3, 4, 8, 10)]
    alphas = [0.9, 0.35, 1.0, 0.12, 0.6, -0.8, -0.25, -1.0, -0.4, -0.72]
    files["sv_svseq.fa"] = fasta(sv)
    files["sv_svalpha.out"] = "".join(f"{n}\t{a:11.6e}\n" for (n, _), a in zip(sv, alphas))
    # --- added in Phase 1 (new PRNG draws come after everything above, so older fixtures are unchanged)
    # the same SV model with 120-character names, and long-named test sequences (name-buffer overflow)
    ln = lambda n: "N" * 110 + "_" + n
    files["sv_longnames_svseq.fa"] = fasta([(ln(n), s) for n, s in sv])
    files["sv_longnames_svalpha.out"] = "".join(f"{ln(n)}\t{a:11.6e}\n" for (n, _), a in zip(sv, alphas))
    files["test_longnames.fa"] = fasta([("T" * 118 + f"_{i}", s) for i, (_, s) in enumerate(
        [(f"t{i}", s) for i, s in enumerate(seqs(rng, 10, 100, motif=MOTIF))])])
    # a fixed SV model over the two-letter alphabet (classify with -A; exercises the b != 4 weights)
    b2p = [(f"b2p{i}", s) for i, s in enumerate(seqs(rng, 10, 80, alpha="AB", motif="ABBABAA"))]
    b2n = [(f"b2n{i}", s) for i, s in enumerate(seqs(rng, 10, 80, alpha="AB"))]
    svb2 = [b2p[i] for i in (0, 3, 6, 9)] + [b2n[i] for i in (1, 4, 7)]
    ab2 = [0.7, 1.0, 0.2, 0.55, -0.9, -0.3, -1.0]
    files["sv_b2_svseq.fa"] = fasta(svb2)
    files["sv_b2_svalpha.out"] = "".join(f"{n}\t{a:11.6e}\n" for (n, _), a in zip(svb2, ab2))
    files["test_b2.fa"] = fasta([(f"tb{i}", s) for i, s in enumerate(seqs(rng, 8, 80, alpha="AB", motif="ABBABAA"))])
    # --- added in Phase 3: name edge cases (kernel values depend only on the sequences)
    edge = [("", small_pos[0][1]),                       # empty name
            ("a\tb c", small_pos[1][1]),                  # tab and space in the header: name is "a"
            ("s\u00e9q_\u00fc\u00f1", small_pos[2][1]),  # non-ASCII (UTF-8) name
            ("X" * 200, small_pos[3][1]),                 # 200-character name
            ("neg0", small_pos[4][1]),                    # same name as a record in the negative file
            ("neg0", small_pos[5][1]),                    # ... twice
            ("pos6", small_pos[6][1]), ("pos7", small_pos[7][1])]
    files["names_edge_pos.fa"] = fasta(edge)
    # SV models that must be rejected: a name listed twice, a name absent from the sequence file
    files["sv_dupalpha_svalpha.out"] = files["sv_svalpha.out"] + "pos2\t 1.000000e-01\n"
    files["sv_missing_svalpha.out"] = files["sv_svalpha.out"] + "pos99\t 1.000000e-01\n"
    # --- added in Phase 5: protein (b=20) sequences, an RNA alphabet file with complement pairs
    AA = "ACDEFGHIKLMNPQRSTVWY"
    files["protein_pos.fa"] = fasta([(f"prp{i}", s) for i, s in enumerate(seqs(rng, 6, 60, alpha=AA, motif="WKRHQ"))])
    files["protein_neg.fa"] = fasta([(f"prn{i}", s) for i, s in enumerate(seqs(rng, 6, 60, alpha=AA))])
    files["alphabet_protein.txt"] = "\n".join(AA) + "\n"
    files["alphabet_rna_pairs.txt"] = "# RNA with complement pairs\nA U\nC G\nG C\nU A\n"
    files["rna_pos.fa"] = fasta([(n, s.replace("T", "U")) for n, s in small_pos[:6]])
    files["rna_neg.fa"] = fasta([(n, s.replace("T", "U")) for n, s in small_neg[:6]])
    files["alphabet_b5_pairs_bad.txt"] = "A T\nC G\nG C\nT A\nN\n"
    files["sv_protein_svalpha.out"] = "".join(f"prn{i}\t{a:11.6e}\n" for i, a in enumerate([0.5, -0.3, 0.8, -1.0, 0.2, -0.6]))
    return files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="verify the committed fixtures match; do not write")
    a = ap.parse_args()
    files = build()
    os.makedirs(FIX, exist_ok=True)
    bad = 0
    for name, content in sorted(files.items()):
        path = os.path.join(FIX, name)
        if a.check:
            cur = open(path, "rb").read() if os.path.exists(path) else None
            if cur != content.encode():
                print(f"MISMATCH {path}")
                bad += 1
        else:
            with open(path, "wb") as f:
                f.write(content.encode())
    if a.check:
        print("fixtures OK" if not bad else f"{bad} fixture(s) differ")
        sys.exit(1 if bad else 0)
    print(f"wrote {len(files)} fixtures to {FIX}")


if __name__ == "__main__":
    main()
