# Engreitz WTC11 Endothelial

**Lab**: Engreitz  •  **Cell line**: WTC11  •  **Differentiation**: Endothelial

Production dataset for the [2026 UTSW jamboree](../../README.md). Currently a placeholder — the data is not on the IGVF portal yet.

## Identity

| | |
|---|---|
| Dataset ID | `Engreitz_WTC11-endothelial-cells_TF-Perturb-seq` |
| IGVF analysis set | _none yet_ |
| Construct library | _none yet_ |
| Guide file | _none yet_ |
| Guide pools | _?_ |
| Perturbation | CRISPRi |
| Assay | CC Perturb-seq (Engreitz lab; same technology as the Engreitz_WTC11-benchmark dataset) |
| Multiplexing | _?_ |
| Measurement sets | 0 |

## Status

| Output | Status | Synapse |
|---|---|---|
| CRISPR pipeline | ☐ blocked: no portal data | — |
| cNMF | ☐ blocked: no upstream data | — |
| Energy distance | ☐ blocked: no upstream data | — |

Open issue: [Engreitz endothelial data not on portal](../../issues/engreitz-no-data.md).

## Source data

None yet. Once the Engreitz team uploads to the IGVF portal under the `TF Perturb-seq Project` collection, the [`scripts/query_igvf_portal.py`](../../scripts/query_igvf_portal.py) snapshot will pick it up and the rest of the pipeline can begin.

## Notes

- Listed in [`2026_05_07_state.png`](../../2026_05_07_state.png) but no IGVF portal entry as of the latest portal snapshot.
- No entry under the local `tf_perturb_seq/datasets/` directory yet either.
- The Engreitz WTC11 benchmark dataset (separate, uses CC Perturb-seq) is the closest reference for what the production pipeline parameters will look like.
