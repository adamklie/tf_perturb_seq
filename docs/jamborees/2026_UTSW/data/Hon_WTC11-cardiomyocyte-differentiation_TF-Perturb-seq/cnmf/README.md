# Hon WTC11 Cardiomyocyte — cNMF

**Status**: ⏳ Not on jamboree Synapse yet. Two efforts in flight:

| Effort | Owner | Source MuData | Status |
|---|---|---|---|
| `051126_honcm_torchcnmf_KskillA` | Adam | `seqspec_v3` (our GCS run) | ❌ Stage 1 Convert script **failed** with a numba JIT working-dir bug (see SLURM err `seqspec_v3/cnmf/Data/convert_10731567.err`). Held until fix or until Alexandra's run lands. |
| (in parallel) | Alexandra | Weizhou's MuData (Synapse [`syn74522725`](https://www.synapse.org/Synapse:syn74522725)) | 🔄 Running on her side. Canonical target for jamboree — cNMF should be on Weizhou's data to stay consistent with the ED runs. |

Tracking: [Issue 20](https://github.com/adamklie/tf_perturb_seq/issues/20).

## Plan once Alexandra delivers

Mirror her cNMF output to `2026_UTSW/datasets/Hon_WTC11-cardiomyocyte-differentiation_TF-Perturb-seq/cnmf/<run_name>/` per the [`schemas/cnmf.json`](../../../schemas/cnmf.json) bundle inclusion rule. Use [`data/mirror_cnmf_outputs.py`](../../../data/mirror_cnmf_outputs.py).

## Schema + walkthrough

- Machine-readable schema: [`schemas/cnmf.json`](../../../schemas/cnmf.json)
- Analysis-level walkthrough: [`docs/analysis/cNMF_OUTPUTS.md`](../../../../../analysis/cNMF_OUTPUTS.md), [`docs/analysis/cNMF.md`](../../../../../analysis/cNMF.md)
