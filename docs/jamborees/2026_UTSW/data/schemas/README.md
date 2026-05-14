# Schemas

JSON schemas for every output table mirrored to the jamboree, in a Frictionless-Table-Schema-inspired format extended with `source` / `note` / `produced_by` keys.

| Schema | Output | Rows | Cols |
|---|---|---|---|
| `tf_metadata.json` | Comprehensive TF metadata | 1,983 | 16 |
| `tf_metadata_simplified.json` | Simplified TF metadata | 1,983 | 8 |
| `experimental_metadata.json` | Comprehensive experimental metadata | 5 | 26 |
| `experimental_metadata_simplified.json` | Simplified experimental metadata | 5 | 12 |
| `guide_metadata.json` | IGVF guide library file (pools A-D), used as-is from the portal | 14,150 | 18 |
| `crispr_pipeline.json` | Per-dataset CRISPR pipeline outputs (3 dirs mirrored as-is) | per-dataset | per-dataset |
| `energy_distance.json` | Per-dataset energy-distance pipeline outputs | per-dataset | per-dataset |
| `cnmf.json` | Per-dataset cNMF run outputs (selected-k plus sweep-as-provenance) | per-dataset | per-dataset |

## Conventions

- **Field types**: `string`, `integer`, `boolean` (Frictionless types).
- **`source`** (custom): the upstream source of a column (a file path, URL, or pipeline-config key).
- **Marker conventions** in tabular outputs are documented per-schema; common values:
  - `?` = value not yet known.
  - `-` = not applicable.
  - `0` (in numeric columns) = unknown.
- **Comprehensive vs simplified**: `tf_metadata` and `experimental_metadata` each have a comprehensive (machine-readable) form and a simplified (human-readable) form. `guide_metadata` is published as-is from the IGVF portal — no simplified version.

## Adding a new schema

1. Create `schemas/<name>.json` describing the output.
2. Reference `produced_by` (the script that builds the output) and list `sources`.
3. Update this README's table.
4. Update the top-level README's outputs section.
