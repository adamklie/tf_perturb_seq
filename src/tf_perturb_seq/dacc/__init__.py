"""DACC (Data and Administrative Coordinating Center) file generators.

Produces the file types defined in the IGVF CPN FG file-format spec from our
canonical inputs (guide library, cNMF outputs, perturbo trans-DE). See
`docs/data/DACC.md` for the audit of what's defined and what's submitted.

Modules:
    build_tf_universe — TF Universe TSV from a harmonized guide library
    build_element_universe — Element Universe BED from a harmonized guide library
"""
