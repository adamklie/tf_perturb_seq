#!/usr/bin/env python3
"""
filter_mudata_by_sgrna_count.py

Removes cells with more than a given number of total assigned sgRNAs
from a MuData object.

Usage:
    python filter_mudata_by_sgrna_count.py \
        --input path/to/input.h5mu \
        --output path/to/output.h5mu \
        [--max-guides 15] \
        [--log path/to/output.log]

    srun \
        --job-name=filter_sgrna \
        --mem=300G \
        --cpus-per-task=4 \
        --time=02:00:00 \
        --partition=common \
        python filter_mudata_by_sgrna_count.py \
            --input /hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/inference_mudata.h5mu \
            --output /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata.h5mu \
            --max-guides 15 \
            --log /hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/logs/filter_sgrna.log
"""

import argparse
import gc
import logging
import sys
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import mudata as md


def setup_logger(log_path: Optional[str]) -> logging.Logger:
    logger = logging.getLogger("filter_sgrna")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if log_path:
        fh = logging.FileHandler(log_path)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def get_guide_counts_h5py(h5mu_path: str, logger: logging.Logger) -> np.ndarray:
    """
    Read ONLY the guide_assignment layer from the h5mu file via h5py.
    Handles CSR, CSC, and dense encodings.
    Returns a 1-D int32 array of assigned-guide counts per cell.
    Peak memory: O(nnz) for CSR indptr only — far below loading full MuData.
    """
    layer_path = "mod/guide/layers/guide_assignment"

    with h5py.File(h5mu_path, "r") as f:
        if layer_path not in f:
            raise KeyError(
                f"Expected layer at '{layer_path}' not found in {h5mu_path}. "
                f"Top-level keys: {list(f.keys())}"
            )
        layer = f[layer_path]

        # Decode encoding-type attribute (bytes in older h5ad versions)
        enc = layer.attrs.get("encoding-type", "")
        if isinstance(enc, bytes):
            enc = enc.decode()
        enc = enc.lower()
        logger.info(f"guide_assignment encoding: '{enc}'")

        if "csr" in enc:
            # Row sums of a binary CSR = diff of indptr
            indptr = layer["indptr"][:]
            counts = np.diff(indptr).astype(np.int32)

        elif "csc" in enc:
            # Each entry in 'indices' is a row index; bincount gives row sums
            shape = tuple(
                layer.attrs["shape"] if "shape" in layer.attrs else layer["shape"][:]
            )
            indices = layer["indices"][:]
            counts = np.bincount(indices, minlength=shape[0]).astype(np.int32)

        else:
            # Dense fallback — read in row-chunks to limit peak memory
            logger.warning(
                f"Unrecognised encoding '{enc}'; falling back to chunked dense read."
            )
            data = layer["data"] if "data" in layer else layer
            n_rows = data.shape[0]
            counts = np.zeros(n_rows, dtype=np.int32)
            chunk = 50_000
            for start in range(0, n_rows, chunk):
                block = data[start : start + chunk]
                counts[start : start + chunk] = (block != 0).sum(axis=1)

    return counts


def build_keep_mask(
    guide_counts: np.ndarray,
    max_guides: int,
    logger: logging.Logger,
) -> np.ndarray:
    n_cells = len(guide_counts)

    assigned_mask = guide_counts > 0
    n_with_guides = int(assigned_mask.sum())
    logger.info(
        f"Cells with at least 1 assigned guide: "
        f"{n_with_guides:,} / {n_cells:,} "
        f"({100 * n_with_guides / n_cells:.1f}%)"
    )
    logger.info(
        f"Assigned guides per cell (all cells):  "
        f"median={int(np.median(guide_counts))}, "
        f"mean={guide_counts.mean():.2f}, "
        f"max={int(guide_counts.max())}"
    )

    keep_mask = guide_counts <= max_guides
    n_removed = int((~keep_mask).sum())
    n_kept = int(keep_mask.sum())

    logger.info(f"Threshold: max_guides = {max_guides}")
    logger.info(
        f"Cells removed  (> {max_guides} guides): "
        f"{n_removed:,} / {n_cells:,} ({100 * n_removed / n_cells:.2f}%)"
    )
    logger.info(
        f"Cells retained (<= {max_guides} guides): "
        f"{n_kept:,} / {n_cells:,} ({100 * n_kept / n_cells:.2f}%)"
    )

    return keep_mask


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter MuData cells by total assigned sgRNA count."
    )
    parser.add_argument("--input",  "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--max-guides", "-m", type=int, default=15)
    parser.add_argument("--log", "-l", default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logger = setup_logger(args.log)

    input_path  = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Step 1 — compute mask cheaply (h5py only, no MuData load yet)
    # ------------------------------------------------------------------ #
    logger.info("Computing guide counts via h5py (low-memory pass)...")
    guide_counts = get_guide_counts_h5py(str(input_path), logger)
    keep_mask    = build_keep_mask(guide_counts, args.max_guides, logger)
    del guide_counts
    gc.collect()

    # ------------------------------------------------------------------ #
    # Step 2 — load full MuData (one copy only)
    # ------------------------------------------------------------------ #
    logger.info(f"Reading full MuData from: {input_path}")
    mdata = md.read_h5mu(str(input_path))
    logger.info(f"Loaded: {mdata.n_obs:,} cells × {mdata.n_vars:,} vars")

    # ------------------------------------------------------------------ #
    # Step 3 — filter as a VIEW and write directly (no .copy())
    # ------------------------------------------------------------------ #
    logger.info("Filtering (view, no copy)...")
    mdata_view = mdata[keep_mask]

    # Free the unfiltered data before writing to avoid 2× peak memory
    del mdata
    gc.collect()

    logger.info(f"Writing filtered MuData to: {output_path}")
    mdata_view.write_h5mu(str(output_path))

    # Post-filter summary
    logger.info(
        f"Output MuData: {mdata_view.n_obs:,} cells × {mdata_view.n_vars:,} vars"
    )
    logger.info("Done.")


if __name__ == "__main__":
    main()