#!/usr/bin/env python3
"""Filter the GTDB SSU reference by minimum length, then take a random
holdout (default 10 %) as a held-out test set for BLCA benchmarking.

The output layout is what the BLCA Nextflow taxonomy_only entry needs:
  reference_90.fasta       (90 % retained sequences — train reference)
  reference_90.taxonomy    (matching accession → lineage map)
  test_10.fasta            (10 % held-out queries)
  test_10.taxonomy         (ground truth lineages for the held-out set)
  split_stats.json         (counts of the split)

We deliberately do NOT stratify by species: a pure random split, as manuscripts
usually report. Test species absent from the training fraction are scored one rank
up by the scoring scripts (analysis/01_holdout_accuracy.R).

Usage:
  python split_gtdb_holdout.py \\
      --in-fasta    db/gtdb_ssu_BLCAparsed.fasta \\
      --in-taxonomy db/gtdb_ssu_BLCAparsed.taxonomy \\
      --min-length  1000 \\
      --holdout-frac 0.10 \\
      --seed         42 \\
      --out-dir      benchmark_blca_holdout
"""
import argparse
import json
import random
from pathlib import Path




def read_fasta_records(path):
    """Yield (id, seq_lines) pairs from a fasta file, streaming.

    Avoids loading the whole 1.4 GB file twice. Sequences are returned as
    a list of stripped lines so the caller can `''.join(lines)` only when
    it needs to count length or write — never doubling RAM for unkept seqs.
    """
    rec_id = None
    rec_lines = []
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
    ap.add_argument("--min-length",   type=int,   default=1000)
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

    # PASS 1: find all accession IDs whose sequence is >= min-length.
    # Streaming, so RAM stays bounded.
    print(f"[split] pass 1: collect accessions with length >= {args.min_length}")
    kept_ids = []
    n_total = 0
    for acc, seq_lines in read_fasta_records(args.in_fasta):
        n_total += 1
        seq_len = sum(len(s) for s in seq_lines)
        if seq_len >= args.min_length:
            kept_ids.append(acc)
    print(f"[split]   {n_total:,} total sequences in input fasta")
    print(f"[split]   {len(kept_ids):,} retained at length >= {args.min_length}  "
          f"({len(kept_ids)/n_total*100:.1f}% kept)")

    # Random 10/90 split (seeded).
    rng = random.Random(args.seed)
    rng.shuffle(kept_ids)
    n_test = int(args.holdout_frac * len(kept_ids))
    test_set = set(kept_ids[:n_test])
    train_set = set(kept_ids[n_test:])
    print(f"[split]   test (10%) : {n_test:,}")
    print(f"[split]   train (90%): {len(train_set):,}")

    # PASS 2: stream-write TRAIN; buffer TEST, then write TEST in the shuffled order of
    # kept_ids[:n_test] rather than the input's taxonomy-grouped order (this is the
    # record order of the published test_10 files). Buffering 10% of GTDB takes ~140 MB.
    test_fa   = out_dir / "test_10.fasta"
    train_fa  = out_dir / "reference_90.fasta"
    test_txt  = out_dir / "test_10.taxonomy"
    train_txt = out_dir / "reference_90.taxonomy"

    print("[split] pass 2: write reference_90 (streamed) + test_10 (shuffled order)")
    n_train_written = n_train_tax_written = 0
    test_seqs = {}
    with open(train_fa, "w") as rf, open(train_txt, "w") as rtx:
        for acc, seq_lines in read_fasta_records(args.in_fasta):
            if acc in test_set:
                test_seqs[acc] = seq_lines
            elif acc in train_set:
                rf.write(f">{acc}\n")
                for s in seq_lines: rf.write(s + "\n")
                n_train_written += 1
                if acc in tax:
                    rtx.write(f"{acc}\t{tax[acc]}\n")
                    n_train_tax_written += 1
            # else: filtered out by length
    n_test_written = n_test_tax_written = 0
    with open(test_fa, "w") as tf, open(test_txt, "w") as ttx:
        for acc in kept_ids[:n_test]:              # shuffled order
            seq_lines = test_seqs.get(acc)
            if seq_lines is None:
                continue
            tf.write(f">{acc}\n")
            for s in seq_lines: tf.write(s + "\n")
            n_test_written += 1
            if acc in tax:
                ttx.write(f"{acc}\t{tax[acc]}\n")
                n_test_tax_written += 1

    stats = {
        "input": {
            "in_fasta":              args.in_fasta,
            "in_taxonomy":           args.in_taxonomy,
            "n_total_sequences":     n_total,
            "min_length":            args.min_length,
            "n_kept_after_length":   len(kept_ids),
        },
        "split": {
            "seed":                  args.seed,
            "holdout_frac":          args.holdout_frac,
            "n_train":               n_train_written,
            "n_test":                n_test_written,
            "n_train_tax_lines":     n_train_tax_written,
            "n_test_tax_lines":      n_test_tax_written,
        },
    }
    with open(out_dir / "split_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)

    print("\n" + json.dumps(stats, indent=2))
    print(f"\n[split] wrote outputs to {out_dir}")


if __name__ == "__main__":
    main()
