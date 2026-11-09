#!/usr/bin/env python3
"""Random 90/10 holdout split of the BLCA-parsed UNITE reference for the fungal
ITS classifier benchmark. Fungal analogue of split_gtdb_holdout.py.

Identical logic to the GTDB version (a pure, seeded random split, NOT stratified
by species).
Two differences only:
  * default --min-length is 0 (ITS is short; UNITE median ~526 bp). The 1000 bp
    16S threshold would discard the entire database. Length filtering, if any,
    is normally done once in build_unite_blca_db.sh.
  * messages/labels say UNITE/ITS.

Outputs (into --out-dir):
  reference_90.fasta / reference_90.taxonomy   90% train reference
  test_10.fasta      / test_10.taxonomy        10% held-out queries + ground truth
  split_stats.json                             counts of the split
"""
import argparse
import json
import random
from pathlib import Path


def read_fasta_records(path):
    rec_id, rec_lines = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, rec_lines
                rec_id = line[1:].strip().split()[0]
                rec_lines = []
            else:
                rec_lines.append(line.strip())
        if rec_id is not None:
            yield rec_id, rec_lines



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-fasta",     required=True)
    ap.add_argument("--in-taxonomy",  required=True)
    ap.add_argument("--min-length",   type=int,   default=0,
                    help="drop refs shorter than this (default 0 = keep all; ITS is short)")
    ap.add_argument("--holdout-frac", type=float, default=0.10)
    ap.add_argument("--seed",         type=int,   default=42)
    ap.add_argument("--out-dir",      required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[split] loading taxonomy from {args.in_taxonomy}")
    tax = {}
    with open(args.in_taxonomy) as fh:
        for line in fh:
            if "\t" not in line:
                continue
            acc, lineage = line.rstrip("\n").split("\t", 1)
            tax[acc] = lineage
    print(f"[split]   {len(tax):,} taxonomy entries loaded")

    print(f"[split] pass 1: collect accessions with length >= {args.min_length}")
    kept_ids = []
    n_total = 0
    for acc, seq_lines in read_fasta_records(args.in_fasta):
        n_total += 1
        if sum(len(s) for s in seq_lines) >= args.min_length:
            kept_ids.append(acc)
    print(f"[split]   {n_total:,} total sequences in input fasta")
    print(f"[split]   {len(kept_ids):,} retained at length >= {args.min_length}  "
          f"({len(kept_ids)/max(1,n_total)*100:.1f}% kept)")

    rng = random.Random(args.seed)
    rng.shuffle(kept_ids)
    n_test = int(args.holdout_frac * len(kept_ids))
    test_set = set(kept_ids[:n_test])
    train_set = set(kept_ids[n_test:])
    print(f"[split]   test ({args.holdout_frac*100:.0f}%) : {n_test:,}")
    print(f"[split]   train          : {len(train_set):,}")

    test_fa   = out_dir / "test_10.fasta"
    train_fa  = out_dir / "reference_90.fasta"
    test_txt  = out_dir / "test_10.taxonomy"
    train_txt = out_dir / "reference_90.taxonomy"

    # PASS 2: stream-write TRAIN; buffer TEST, then write TEST in the shuffled order of
    # kept_ids[:n_test] (the record order of the published test_10 files).
    print("[split] pass 2: write reference_90 (streamed) + test_10 (shuffled order)")
    n_train_written = n_train_tax_written = 0
    test_seqs = {}
    with open(train_fa, "w") as rf, open(train_txt, "w") as rtx:
        for acc, seq_lines in read_fasta_records(args.in_fasta):
            if acc in test_set:
                test_seqs[acc] = seq_lines
            elif acc in train_set:
                rf.write(f">{acc}\n")
                for s in seq_lines:
                    rf.write(s + "\n")
                n_train_written += 1
                if acc in tax:
                    rtx.write(f"{acc}\t{tax[acc]}\n")
                    n_train_tax_written += 1
    n_test_written = n_test_tax_written = 0
    with open(test_fa, "w") as tf, open(test_txt, "w") as ttx:
        for acc in kept_ids[:n_test]:              # shuffled order
            seq_lines = test_seqs.get(acc)
            if seq_lines is None:
                continue
            tf.write(f">{acc}\n")
            for s in seq_lines:
                tf.write(s + "\n")
            n_test_written += 1
            if acc in tax:
                ttx.write(f"{acc}\t{tax[acc]}\n")
                n_test_tax_written += 1

    stats = {
        "input": {
            "in_fasta":            args.in_fasta,
            "in_taxonomy":         args.in_taxonomy,
            "n_total_sequences":   n_total,
            "min_length":          args.min_length,
            "n_kept_after_length": len(kept_ids),
        },
        "split": {
            "seed":              args.seed,
            "holdout_frac":      args.holdout_frac,
            "n_train":           n_train_written,
            "n_test":            n_test_written,
            "n_train_tax_lines": n_train_tax_written,
            "n_test_tax_lines":  n_test_tax_written,
        },
    }
    with open(out_dir / "split_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)

    print("\n" + json.dumps(stats, indent=2))
    print(f"\n[split] wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
