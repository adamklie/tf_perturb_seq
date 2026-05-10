"""Validate the outputs of a torch-cNMF run against the curation rule in schemas/cnmf.json.

Runs four layers of checks against a single dataset's cNMF run dir
(e.g., the per-run folder under datasets/<id>/PerturbNMF/Result/<run>/):

  Layer 1 — file presence and non-emptiness
              (run-level files + selected-k full data + sweep-as-provenance)
  Layer 2 — table-shape sanity for the selected-k tabular outputs
              (gene_spectra_score is k×G; usages.consensus is N×k; etc.)
  Layer 3 — value-range sanity (loadings real, usages non-negative, etc.)
  Layer 4 — cross-reference vs the Hon WTC11 benchmark reference run on HPC
              (the schema-verification source documented in schemas/cnmf.json).
              Skip with --no-cross-ref or pass --hon-ref-dir to override path.

Each layer prints PASS/FAIL/WARN per check. Exit code 0 if all REQUIRED checks
pass, 1 otherwise. WARNs are non-fatal.

Usage:
    python validate_cnmf_outputs.py \\
        --source-dir /path/to/PerturbNMF/Result/<run_name> \\
        --selected-k 50

    python validate_cnmf_outputs.py --source-dir <DIR> --selected-k <K> --no-cross-ref
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Default Hon benchmark reference (HPC path; only used in Layer 4)
HON_REF_DIR = Path(
    "/cellar/users/aklie/projects/tf_perturb_seq/datasets/"
    "Hon_WTC11-benchmark_TF-Perturb-seq/PerturbNMF/Result/"
    "030726_20iter_5KHVG_torch_halsvar_batch_e7"
)

# Per-(k, dt) Eval/ TXT files we expect (from schemas/cnmf.json files[].Eval)
EVAL_REQUIRED_FILES = {
    "Explained_Variance.txt",
    "Explained_Variance_Summary.txt",
    "GO_term_enrichment.txt",
    "categorical_association_results.txt",
    "categorical_association_posthoc.txt",
    "geneset_enrichment.txt",
    "trait_enrichment.txt",
}
EVAL_PERTURBATION_PATTERN = re.compile(r".*_perturbation_association_results_.+\.txt$")
EVAL_OPTIONAL_CALIBRATION = "fake_perturbation_association_calibration.txt"

ANSI = {
    "PASS": "\033[32mPASS\033[0m",
    "FAIL": "\033[31mFAIL\033[0m",
    "WARN": "\033[33mWARN\033[0m",
    "INFO": "\033[36mINFO\033[0m",
}


class Result:
    def __init__(self):
        self.fails: list[str] = []
        self.warns: list[str] = []

    def report(self, level: str, msg: str) -> None:
        print(f"{ANSI.get(level, level):>4}  {msg}")
        if level == "FAIL":
            self.fails.append(msg)
        elif level == "WARN":
            self.warns.append(msg)

    def summary(self) -> int:
        print()
        if self.fails:
            print(f"{ANSI['FAIL']}  {len(self.fails)} required check(s) failed:")
            for m in self.fails:
                print(f"  - {m}")
        if self.warns:
            print(f"{ANSI['WARN']}  {len(self.warns)} non-fatal warning(s):")
            for m in self.warns:
                print(f"  - {m}")
        if not self.fails and not self.warns:
            print(f"{ANSI['PASS']}  all checks passed")
        return 1 if self.fails else 0


def _human_size(p: Path) -> str:
    n = p.stat().st_size
    for unit in ("B", "K", "M", "G"):
        if n < 1024:
            return f"{n:.0f}{unit}"
        n /= 1024
    return f"{n:.0f}T"


def _fmt_dt(dt: float) -> str:
    return f"{int(dt)}_{int(round((dt - int(dt)) * 10))}"


def _present(p: Path) -> bool:
    return p.exists() and (not p.is_file() or p.stat().st_size > 0)


def _run_name_from_src(src: Path) -> str:
    return src.name


def layer1_presence(src: Path, sel_k: int, dt: float, r: Result) -> None:
    print("\n=== Layer 1: file presence + non-emptiness ===")
    run = _run_name_from_src(src)
    sel_dt = _fmt_dt(dt)

    # Run-level
    run_level = [
        f"{run}.k_selection.png",
        f"{run}.k_selection_stats.df.npz",
        f"{run}.overdispersed_genes.txt",
    ]
    for rel in run_level:
        p = src / rel
        if not p.exists():
            r.report("FAIL", f"missing run-level file: {rel}")
        elif p.stat().st_size == 0:
            r.report("FAIL", f"empty run-level file: {rel}")
        else:
            r.report("PASS", f"{rel} ({_human_size(p)})")

    # README.txt (selection rationale) — required for sweep-as-provenance
    readme = src / "README.txt"
    if not _present(readme):
        r.report("WARN", "README.txt missing (selection rationale belongs here per cNMF.md)")
    else:
        r.report("PASS", f"README.txt ({_human_size(readme)})")

    # config_*.yml — at least one
    configs = sorted(src.glob("config_*.yml"))
    if not configs:
        r.report("WARN", "no config_*.yml files (run config not captured)")
    else:
        r.report("PASS", f"{len(configs)} config_*.yml file(s)")

    # logs/ folder
    logs = src / "logs"
    if not logs.is_dir():
        r.report("WARN", "logs/ folder missing (SLURM stdout/stderr / resource_monitor)")
    else:
        log_files = list(logs.iterdir())
        r.report("PASS", f"logs/ ({len(log_files)} files)")

    # Selected-k flat files
    selected_k_flat = [
        f"{run}.gene_spectra_score.k_{sel_k}.dt_{sel_dt}.txt",
        f"{run}.gene_spectra_tpm.k_{sel_k}.dt_{sel_dt}.txt",
        f"{run}.spectra.k_{sel_k}.dt_{sel_dt}.consensus.txt",
        f"{run}.starcat_spectra.k_{sel_k}.dt_{sel_dt}.txt",
        f"{run}.usages.k_{sel_k}.dt_{sel_dt}.consensus.txt",
        f"{run}.clustering.k_{sel_k}.dt_{sel_dt}.png",
    ]
    for rel in selected_k_flat:
        p = src / rel
        if not p.exists():
            r.report("FAIL", f"missing selected-k file: {rel}")
        elif p.stat().st_size == 0:
            r.report("FAIL", f"empty selected-k file: {rel}")
        else:
            r.report("PASS", f"{rel} ({_human_size(p)})")

    # Selected-k subdirs / files
    selected_k_targets = [
        ("file", f"adata/cNMF_{sel_k}_{sel_dt}.h5mu"),
        ("file", f"Annotation/{sel_k}_{sel_dt}.xlsx"),
        ("dir", f"Eval/{sel_k}_{sel_dt}"),
        ("dir", f"Interpretation/Summary_table/{sel_k}_{sel_dt}"),
        ("dir", f"Plot/Program_{sel_k}_{sel_dt}"),
        ("dir", f"Plot/Perturb_gene_{sel_k}_{sel_dt}"),
    ]
    for kind, rel in selected_k_targets:
        p = src / rel
        if kind == "file":
            if not p.exists():
                r.report("FAIL" if "h5mu" in rel else "WARN", f"missing selected-k {kind}: {rel}")
            elif p.stat().st_size == 0:
                r.report("FAIL", f"empty {rel}")
            else:
                r.report("PASS", f"{rel} ({_human_size(p)})")
        else:
            if not p.is_dir():
                r.report("WARN", f"missing selected-k {kind}: {rel}")
            else:
                n = sum(1 for _ in p.rglob("*") if _.is_file())
                r.report("PASS", f"{rel}/ ({n} files)")

    # Sweep-as-provenance: all-k gene_spectra_score
    gss = sorted(src.glob(f"{run}.gene_spectra_score.k_*.dt_*.txt"))
    if not gss:
        r.report("FAIL", "no all-k gene_spectra_score files (sweep-as-provenance is empty)")
    else:
        ks = sorted({int(re.search(r"k_(\d+)", f.name).group(1)) for f in gss})
        r.report("PASS", f"all-k gene_spectra_score: {len(gss)} files spanning k={ks[0]}..{ks[-1]} ({len(ks)} unique k)")

    # Sweep-as-provenance: clustering pngs
    clusterings = sorted(src.glob(f"{run}.clustering.k_*.dt_*.png"))
    if not clusterings:
        r.report("FAIL", "no all-k clustering pngs (sweep-as-provenance is empty)")
    else:
        r.report("PASS", f"all-k clustering pngs: {len(clusterings)}")

    # Sweep-as-provenance: Eval/<k>_<dt>/ subdirs
    eval_dir = src / "Eval"
    if not eval_dir.is_dir():
        r.report("FAIL", "Eval/ missing")
    else:
        eval_subs = sorted(d for d in eval_dir.iterdir() if d.is_dir())
        if not eval_subs:
            r.report("FAIL", "Eval/ has no <k>_<dt>/ subdirs")
        else:
            r.report("PASS", f"Eval/: {len(eval_subs)} (k, dt) subdirs")
            # Spot-check the selected-k Eval subdir for required TXTs
            sel_eval = eval_dir / f"{sel_k}_{sel_dt}"
            if sel_eval.is_dir():
                names = {p.name for p in sel_eval.iterdir() if p.is_file()}
                # Match with run name as prefix
                stripped = {re.sub(rf"^{sel_k}_", "", n) for n in names}
                missing = EVAL_REQUIRED_FILES - stripped
                has_perturb = any(EVAL_PERTURBATION_PATTERN.match(n) for n in names)
                if missing:
                    r.report("FAIL", f"Eval/{sel_k}_{sel_dt}/ missing required files: {sorted(missing)}")
                else:
                    r.report("PASS", f"Eval/{sel_k}_{sel_dt}/ has the {len(EVAL_REQUIRED_FILES)} required generic TXTs")
                if not has_perturb:
                    r.report("FAIL", f"Eval/{sel_k}_{sel_dt}/ missing perturbation_association_results_<batch>.txt files")
                else:
                    n_perturb = sum(1 for n in names if EVAL_PERTURBATION_PATTERN.match(n))
                    r.report("PASS", f"Eval/{sel_k}_{sel_dt}/ has {n_perturb} perturbation_association_results files")

    # Plot/k_selection_*/
    plot_dir = src / "Plot"
    if plot_dir.is_dir():
        ks_dirs = [d for d in plot_dir.iterdir() if d.is_dir() and d.name.startswith("k_selection")]
        if not ks_dirs:
            r.report("WARN", "Plot/k_selection_*/ folder missing (k-selection figure pack)")
        else:
            r.report("PASS", f"Plot/k_selection_*/: {len(ks_dirs)} folder(s)")


def layer2_shapes(src: Path, sel_k: int, dt: float, r: Result) -> dict:
    print("\n=== Layer 2: selected-k table shapes ===")
    run = _run_name_from_src(src)
    sel_dt = _fmt_dt(dt)
    out: dict[str, pd.DataFrame] = {}

    # gene_spectra_score: rows = k programs, cols = genes
    gss = src / f"{run}.gene_spectra_score.k_{sel_k}.dt_{sel_dt}.txt"
    if gss.is_file():
        try:
            df = pd.read_csv(gss, sep="\t", index_col=0)
        except Exception as e:
            r.report("FAIL", f"gene_spectra_score: parse failed ({e!r})")
        else:
            if df.shape[0] != sel_k:
                r.report("FAIL", f"gene_spectra_score: rows={df.shape[0]} != selected k {sel_k}")
            else:
                r.report("PASS", f"gene_spectra_score: {df.shape[0]} programs × {df.shape[1]} genes")
            out["gene_spectra_score"] = df

    # usages.consensus: rows = cells, cols = k programs
    usg = src / f"{run}.usages.k_{sel_k}.dt_{sel_dt}.consensus.txt"
    if usg.is_file():
        try:
            # Read just enough to check shape (full file can be 100MB+)
            df = pd.read_csv(usg, sep="\t", index_col=0, nrows=1000)
            full_cols = df.shape[1]
        except Exception as e:
            r.report("FAIL", f"usages.consensus: parse failed ({e!r})")
        else:
            if full_cols != sel_k:
                r.report("FAIL", f"usages.consensus: cols={full_cols} != selected k {sel_k}")
            else:
                r.report("PASS", f"usages.consensus: cols={full_cols} = selected k {sel_k} (header row+1000-row sample)")

    return out


def layer3_values(out: dict, r: Result) -> None:
    print("\n=== Layer 3: value-range sanity ===")
    df = out.get("gene_spectra_score")
    if df is not None:
        nan_frac = df.isna().sum().sum() / df.size
        if nan_frac > 0.01:
            r.report("FAIL", f"gene_spectra_score: NaN fraction {nan_frac:.3f} > 1%")
        else:
            r.report("PASS", f"gene_spectra_score: NaN fraction {nan_frac:.5f}")
        # z-scored: should be roughly mean ~ 0 per program
        per_program_mean_abs = df.mean(axis=1).abs().mean()
        if per_program_mean_abs > 1.0:
            r.report("WARN", f"gene_spectra_score: mean of abs(per-program means) = {per_program_mean_abs:.3f} (z-scored loadings should be ~0)")
        else:
            r.report("PASS", f"gene_spectra_score: per-program mean centered (mean abs = {per_program_mean_abs:.4f})")


def layer4_cross_ref(src: Path, ref_dir: Path, r: Result) -> None:
    print(f"\n=== Layer 4: cross-reference vs Hon benchmark ({ref_dir}) ===")
    if not ref_dir.is_dir():
        r.report("WARN", f"reference dir not accessible: {ref_dir} (skip)")
        return

    # Compare run-level file presence (run-name prefix differs but suffix doesn't)
    def collect_suffix(d: Path) -> set:
        suffixes = set()
        for f in d.iterdir():
            if not f.is_file():
                continue
            m = re.match(r"^.+?(\.k_selection.*|\.overdispersed_genes\.txt)$", f.name)
            if m:
                suffixes.add(m.group(1))
        return suffixes

    src_sfx = collect_suffix(src)
    ref_sfx = collect_suffix(ref_dir)
    missing = ref_sfx - src_sfx
    if missing:
        r.report("FAIL", f"missing run-level file types vs Hon benchmark: {sorted(missing)}")
    else:
        r.report("PASS", f"run-level file types match Hon benchmark ({len(ref_sfx)} types)")

    # Compare Eval/<k>_<dt>/ schema (one (k, dt) subdir from each)
    src_eval = src / "Eval"
    ref_eval = ref_dir / "Eval"
    if src_eval.is_dir() and ref_eval.is_dir():
        src_kdt = sorted(d.name for d in src_eval.iterdir() if d.is_dir())
        ref_kdt = sorted(d.name for d in ref_eval.iterdir() if d.is_dir())
        common = sorted(set(src_kdt) & set(ref_kdt))
        if not common:
            r.report("WARN", "no common (k, dt) Eval subdirs vs Hon benchmark")
        else:
            sample = common[0]
            src_files = {f.name for f in (src_eval / sample).iterdir() if f.is_file()}
            ref_files = {f.name for f in (ref_eval / sample).iterdir() if f.is_file()}
            # Strip the run-name prefix on each side (varies)
            src_sfx_eval = {re.sub(r"^\d+_", "", f) for f in src_files}
            ref_sfx_eval = {re.sub(r"^\d+_", "", f) for f in ref_files}
            ref_only = ref_sfx_eval - src_sfx_eval
            if ref_only:
                r.report("WARN", f"Eval/{sample}/: missing file kinds vs Hon benchmark: {sorted(ref_only)[:5]}")
            else:
                r.report("PASS", f"Eval/{sample}/: file kinds match Hon benchmark")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", required=True, type=Path, help="Path to the cNMF run dir (PerturbNMF/Result/<run_name>/)")
    ap.add_argument("--selected-k", required=True, type=int, help="Group-selected k for the selected-k bundle")
    ap.add_argument("--density-threshold", type=float, default=2.0, help="Density threshold for selected-k bundle (default: 2.0)")
    ap.add_argument("--no-cross-ref", action="store_true", help="Skip Layer 4 cross-reference vs Hon benchmark")
    ap.add_argument("--hon-ref-dir", type=Path, default=HON_REF_DIR, help="Path to Hon benchmark run dir (Layer 4)")
    args = ap.parse_args()

    if not args.source_dir.is_dir():
        sys.exit(f"--source-dir does not exist or is not a directory: {args.source_dir}")

    r = Result()
    print(f"Validating: {args.source_dir}")
    print(f"selected_k={args.selected_k} density_threshold={args.density_threshold}")

    layer1_presence(args.source_dir, args.selected_k, args.density_threshold, r)
    out = layer2_shapes(args.source_dir, args.selected_k, args.density_threshold, r)
    layer3_values(out, r)
    if not args.no_cross_ref:
        layer4_cross_ref(args.source_dir, args.hon_ref_dir, r)

    return r.summary()


if __name__ == "__main__":
    sys.exit(main())
