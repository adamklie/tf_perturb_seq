# TFP3 Roadmap

**Target**: UTSW Jamboree, May 13-16, 2026

**GitHub Project**: [TFP3](https://github.com/users/adamklie/projects/4)

## Milestones

### 1. Finalize benchmark CRISPR pipeline outputs (#17)
Get clean, uniform pipeline inputs that feed into everything downstream. Rerun Engreitz with samplesheet fix, resolve pipeline version mismatch, upload finalized outputs to Synapse for collaborators.

**Unblocks**: milestones 2, 3, and 4.

### 2. Energy distance on all benchmarks (#13)
Chikara Takeuchi (Hon Lab) runs the E-dist pipeline on all 5 benchmark datasets. We review, document, and compare results across technologies.

**Status**: 1/5 delivered (Gersbach HTv2). Waiting on Chikara for remaining 4.

### 3. cNMF on all benchmarks (#14)
Alexandra Mo (Stanford) runs cNMF-GPU on all 5 benchmark datasets. We review k-selection, program annotations, and cross-dataset program comparisons.

**Status**: In progress. Outputs shared via Globus (manual transfer).

### 4. Full cross-tech comparison (#15)
The synthesis. Once milestones 1-3 are complete, integrate QC + repression efficiency + DEG calling + energy distance + cNMF into a coherent benchmark story. This is Figure 1 material for the paper.

**Depends on**: milestones 1, 2, and 3.

### 5. Production datasets through CRISPR pipeline (#16)
Get >= 2 production-scale datasets through the standardized pipeline. Team effort — Adam unblocks, others execute.

| Dataset | Owner | Status |
|---------|-------|--------|
| Hon cardiomyocyte | Lucas/Weizhou/Adam | Sample metadata generated, waiting on Mpathi |
| Gersbach hepatocyte | Sara Geraghty | Troubleshooting pipeline |
| Huangfu HUES8 endoderm | Gary Yang | Blocked on portal analysis set |
| Huangfu HUES8 stemcell | Gary Yang | Blocked on portal analysis set |
| Engreitz endothelial | Olga (Engreitz Lab) | TBD |

## Dependencies

```
[1] Finalize benchmark pipeline outputs
 |         |         |
 v         v         v
[2] E-dist [3] cNMF [5] Production (independent)
 |         |
 v         v
[4] Full cross-tech comparison → Paper Figure 1
```

Milestones 1-3 can run in parallel. Milestone 4 is the synthesis. Milestone 5 is independent.
