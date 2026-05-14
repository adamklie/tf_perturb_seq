# Reference tables

Cross-dataset reference tables used by working groups and participant analyses. The tables here are computed once across the 5 production datasets and re-used during the jamboree, so per-group analyses don't each re-derive the TF list, gene annotations, or experimental metadata. Bulky raw inputs are kept in this folder too (a GTF and a portal-snapshot CSV) so the generators are runnable without reaching out to other parts of the repo.

Simpler `_simplified` variants of `tf_metadata` and `experimental_metadata` live alongside the full tables — same row keys, fewer columns, suitable for slide-deck cells. Schemas for both variants live in [`../data/schemas/`](../data/schemas/).

## Tables

| File | Description | Generator | Schema |
|---|---|---|---|
| [`tf_metadata.tsv`](tf_metadata.tsv) | One row per target gene in the guide library, with Lambert 2018 DBD / assessment + JASPAR class / family + resolved ENSG | [`scripts/generate_tf_metadata.py`](scripts/generate_tf_metadata.py) | [`tf_metadata.json`](../data/schemas/tf_metadata.json) |
| [`tf_metadata_simplified.tsv`](tf_metadata_simplified.tsv) | Same rows, slide-deck-sized column subset | same | [`tf_metadata_simplified.json`](../data/schemas/tf_metadata_simplified.json) |
| [`experimental_metadata.tsv`](experimental_metadata.tsv) | One row per dataset: lab, cell line, differentiation, IGVF accessions, pipeline run label, parsed config params | [`scripts/generate_experimental_metadata.py`](scripts/generate_experimental_metadata.py) | [`experimental_metadata.json`](../data/schemas/experimental_metadata.json) |
| [`experimental_metadata_simplified.tsv`](experimental_metadata_simplified.tsv) | Same rows, slide-deck-sized column subset | same | [`experimental_metadata_simplified.json`](../data/schemas/experimental_metadata_simplified.json) |
| [`gene_annotations.tsv`](gene_annotations.tsv) | One row per gene from the IGVF GTF: version-stripped `gene_id`, symbol, type, chrom, start, end, strand | [`scripts/build_gene_annotations.py`](scripts/build_gene_annotations.py) | *[FILL IN]* |
| [`gene_disease_associations.tsv`](gene_disease_associations.tsv) | One row per gene, with MONDO + OMIM disease IDs from the HPO `genes_to_disease` release | [`scripts/fetch_hpo_gene_disease.py`](scripts/fetch_hpo_gene_disease.py) | *[FILL IN]* |
| [`tf_universe.tsv`](tf_universe.tsv) | DACC-spec TF universe deliverable for the guide library | *[FILL IN]* | *[FILL IN]* |
| [`element_universe.bed`](element_universe.bed) | DACC-spec element universe deliverable (BED4) for the guide library | *[FILL IN]* | *[FILL IN]* |


## Scripts

All scripts live in [`scripts/`](scripts/) and are run from the jamboree root ([`../`](../)). Outputs land alongside this README.

| Script | Outputs | Example |
|---|---|---|
| [`scripts/generate_tf_metadata.py`](scripts/generate_tf_metadata.py) | `tf_metadata.tsv` + `tf_metadata_simplified.tsv` | `uv run python reference/scripts/generate_tf_metadata.py` |
| [`scripts/generate_experimental_metadata.py`](scripts/generate_experimental_metadata.py) | `experimental_metadata.tsv` + `experimental_metadata_simplified.tsv` | `uv run python reference/scripts/generate_experimental_metadata.py` |
| [`scripts/build_gene_annotations.py`](scripts/build_gene_annotations.py) | `gene_annotations.tsv` | `uv run python reference/scripts/build_gene_annotations.py` |
| [`scripts/fetch_hpo_gene_disease.py`](scripts/fetch_hpo_gene_disease.py) | `gene_disease_associations.tsv` (caches raw download under [`_cache/`](_cache/)) | `uv run python reference/scripts/fetch_hpo_gene_disease.py` |

`generate_tf_metadata.py` and `generate_experimental_metadata.py` read from elsewhere in the repo (`ref/guide_libraries/`, `ref/genome/`, `datasets/<ds>/*.config`); the other two are self-contained against the files in this folder.

## `_cache/`

Raw downloads kept around so generators don't re-fetch on every run (currently only the HPO `genes_to_disease.txt`). Regenerable, gitignored if needed, and not part of the participant-facing surface — don't edit by hand.
