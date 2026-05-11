#\!/bin/bash
set -e

# Patch files for Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3
PATCH_DIR="gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/patch"

echo "Decompressing IGVFFI1697XAXI.tsv"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/07/23/84d7b52d-effd-4227-92ce-0f84bab0e47f/IGVFFI1697XAXI.tsv.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1697XAXI.tsv"

echo "Decompressing IGVFFI5765HMZH.tsv"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/04/cc3d00eb-118c-4230-945f-3cdb847cd849/IGVFFI5765HMZH.tsv.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5765HMZH.tsv"

echo "Decompressing IGVFFI5811SVMQ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/06551b29-c1df-435b-8072-aa8902a88a4d/IGVFFI5811SVMQ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5811SVMQ.yaml"

echo "Decompressing IGVFFI5714WXET.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/08204abc-2674-4beb-ac2a-aad5ae5a3bad/IGVFFI5714WXET.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5714WXET.yaml"

echo "Decompressing IGVFFI8067TGTN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/0a34e56b-d65b-4393-826e-6a661e2f3f5d/IGVFFI8067TGTN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8067TGTN.yaml"

echo "Decompressing IGVFFI9509PDFO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/1239c814-7aa3-4a13-b139-ddacdd5d36fb/IGVFFI9509PDFO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9509PDFO.yaml"

echo "Decompressing IGVFFI4703HIFR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/1806591f-0bc5-45c9-80b8-25f6022b1dc9/IGVFFI4703HIFR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4703HIFR.yaml"

echo "Decompressing IGVFFI7139JADC.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/1844b40e-b37b-442a-b430-de41628f8cd6/IGVFFI7139JADC.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7139JADC.yaml"

echo "Decompressing IGVFFI6820QRAF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/19988181-0b6f-4f0e-b322-3ac7b3eedf9e/IGVFFI6820QRAF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6820QRAF.yaml"

echo "Decompressing IGVFFI1827NALP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/1cd8ec6e-05b5-4a7c-98c7-f877964da21d/IGVFFI1827NALP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1827NALP.yaml"

echo "Decompressing IGVFFI7356TFEN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/22967454-83e2-4ccc-9b7c-25d4e8e2ec49/IGVFFI7356TFEN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7356TFEN.yaml"

echo "Decompressing IGVFFI3077AIXI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/247bb4c6-9851-46d0-a402-db2b7cfcda97/IGVFFI3077AIXI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3077AIXI.yaml"

echo "Decompressing IGVFFI6467WMOP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/4046ca81-a528-4ce2-8706-d65d4fe69bc2/IGVFFI6467WMOP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6467WMOP.yaml"

echo "Decompressing IGVFFI4719TXUG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/43fd6548-bb3b-42e3-ae67-b224938fd164/IGVFFI4719TXUG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4719TXUG.yaml"

echo "Decompressing IGVFFI0762SJJM.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/45eaf625-a6f5-4a05-8dfe-76850350cea6/IGVFFI0762SJJM.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0762SJJM.yaml"

echo "Decompressing IGVFFI0123VGVY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/495a8751-c326-44c6-9f21-155e7a40ab3b/IGVFFI0123VGVY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0123VGVY.yaml"

echo "Decompressing IGVFFI7172DLGO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/4a5206e3-3a6d-493d-871a-378f6a6bd4fa/IGVFFI7172DLGO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7172DLGO.yaml"

echo "Decompressing IGVFFI8609HZFG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/4ad0a93d-475b-4689-a3ff-e45797c0a343/IGVFFI8609HZFG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8609HZFG.yaml"

echo "Decompressing IGVFFI6179UQQW.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/4f22e044-fcc3-441b-8ec3-101bd320bf4c/IGVFFI6179UQQW.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6179UQQW.yaml"

echo "Decompressing IGVFFI9161ALQS.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/4f42b87e-12bf-4a44-af37-8333745cf1ed/IGVFFI9161ALQS.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9161ALQS.yaml"

echo "Decompressing IGVFFI2572SITG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/50f54f49-0cb2-4ba6-82db-7ea91ea0dde3/IGVFFI2572SITG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2572SITG.yaml"

echo "Decompressing IGVFFI3797VTQZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5194af3c-833d-4800-ad50-d6c455a0305b/IGVFFI3797VTQZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3797VTQZ.yaml"

echo "Decompressing IGVFFI1806RVDJ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/528eb764-6334-4008-bde9-b541a3702103/IGVFFI1806RVDJ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1806RVDJ.yaml"

echo "Decompressing IGVFFI4531DYMT.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/54568c2f-ac74-4be7-a90c-71d4c11beb5d/IGVFFI4531DYMT.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4531DYMT.yaml"

echo "Decompressing IGVFFI2498HIKH.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/54f5ae30-f216-44fb-a23a-5708258bc2cd/IGVFFI2498HIKH.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2498HIKH.yaml"

echo "Decompressing IGVFFI0070NUPR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/56652ca9-16ec-472a-9892-0b22215ab805/IGVFFI0070NUPR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0070NUPR.yaml"

echo "Decompressing IGVFFI1362OCGY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/57b9f667-331e-4bce-87b1-1e7cc78ac093/IGVFFI1362OCGY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1362OCGY.yaml"

echo "Decompressing IGVFFI6983XNSI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5acb07b1-8d3d-471c-bbae-22b6612bc48a/IGVFFI6983XNSI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6983XNSI.yaml"

echo "Decompressing IGVFFI5098SJXF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5c90d6c5-ae9c-413e-b0cf-f26989f00b2c/IGVFFI5098SJXF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5098SJXF.yaml"

echo "Decompressing IGVFFI0499BDIV.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5e6519aa-c1a2-4147-b321-f850ebf1c95a/IGVFFI0499BDIV.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0499BDIV.yaml"

echo "Decompressing IGVFFI2726WCUU.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5e7af4e5-dd16-4622-90d3-c5873ac601da/IGVFFI2726WCUU.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2726WCUU.yaml"

echo "Decompressing IGVFFI7302TCCM.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/5ee1a294-8ba5-46c3-b586-727f4bc6253c/IGVFFI7302TCCM.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7302TCCM.yaml"

echo "Decompressing IGVFFI7547UAZY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/639a4009-ee54-4e38-aed7-2dcf2e284321/IGVFFI7547UAZY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7547UAZY.yaml"

echo "Decompressing IGVFFI8405LTTI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/6464b907-ab49-4059-a4d5-e3402ae10825/IGVFFI8405LTTI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8405LTTI.yaml"

echo "Decompressing IGVFFI7752HNLN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/65eef538-0ccc-4a39-8aaa-0e8384d15486/IGVFFI7752HNLN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7752HNLN.yaml"

echo "Decompressing IGVFFI4678MXPS.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/66eb228c-a51d-4f05-86fe-d5c3509d4bad/IGVFFI4678MXPS.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4678MXPS.yaml"

echo "Decompressing IGVFFI2141HZII.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/6969ea5e-cf20-406b-b4e9-75b522b49983/IGVFFI2141HZII.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2141HZII.yaml"

echo "Decompressing IGVFFI7957XWGJ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/6f4045fc-cb7f-4635-9872-b7e8fa4829fe/IGVFFI7957XWGJ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7957XWGJ.yaml"

echo "Decompressing IGVFFI8500HGMG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/6fcc5927-b348-4ffd-bba7-14855a26c75c/IGVFFI8500HGMG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8500HGMG.yaml"

echo "Decompressing IGVFFI9137FWSI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/7322fbc6-b252-4c60-9179-688a41ef8c7d/IGVFFI9137FWSI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9137FWSI.yaml"

echo "Decompressing IGVFFI9685JXEU.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/743e80b2-8ff7-47ac-b778-7e811a9ea9af/IGVFFI9685JXEU.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9685JXEU.yaml"

echo "Decompressing IGVFFI7683SIHY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/76ddf2ce-a7b4-43ec-807f-994c13a92007/IGVFFI7683SIHY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7683SIHY.yaml"

echo "Decompressing IGVFFI4130EZOR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/787cdcc3-f88d-488f-8799-44ab4a62a774/IGVFFI4130EZOR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4130EZOR.yaml"

echo "Decompressing IGVFFI0268RYST.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/790315f9-e21e-43d1-9e86-6cfa5ffa0574/IGVFFI0268RYST.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0268RYST.yaml"

echo "Decompressing IGVFFI7883PGFF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/7aa32b65-8a27-4488-b1ca-29a21e8e7f19/IGVFFI7883PGFF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7883PGFF.yaml"

echo "Decompressing IGVFFI6985XRIB.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/7b40ce57-8b2b-4a81-a320-1d09779768b0/IGVFFI6985XRIB.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6985XRIB.yaml"

echo "Decompressing IGVFFI9977KNWW.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/89e171c0-a859-4384-9e16-2104b60533f2/IGVFFI9977KNWW.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9977KNWW.yaml"

echo "Decompressing IGVFFI7379IGPW.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/92076aad-6f8e-4fbe-984d-dda8d694fb41/IGVFFI7379IGPW.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7379IGPW.yaml"

echo "Decompressing IGVFFI8449BLVI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/97a3ad89-2085-479c-88bb-133d12115765/IGVFFI8449BLVI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8449BLVI.yaml"

echo "Decompressing IGVFFI3001OCFQ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/9894099d-2eac-4def-872c-8d62d954d46f/IGVFFI3001OCFQ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3001OCFQ.yaml"

echo "Decompressing IGVFFI8893TIGY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/9f0587f1-83ba-48c0-bab6-d6cda0f31a28/IGVFFI8893TIGY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8893TIGY.yaml"

echo "Decompressing IGVFFI5990WSFI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/a4557fe5-e096-4762-973b-14acb9d7fa84/IGVFFI5990WSFI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5990WSFI.yaml"

echo "Decompressing IGVFFI5503TDVO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/aaae1c3c-6c17-4cb5-9611-e8ce0e47189f/IGVFFI5503TDVO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5503TDVO.yaml"

echo "Decompressing IGVFFI1335NPXV.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/af5aac0d-a03a-4193-81cf-18f7733c9bdb/IGVFFI1335NPXV.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1335NPXV.yaml"

echo "Decompressing IGVFFI3074VZJW.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/b0738e69-7c67-448e-9b78-0429285509f9/IGVFFI3074VZJW.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3074VZJW.yaml"

echo "Decompressing IGVFFI8791HIOU.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/b198b5a0-7859-4be8-90d3-f0910e6ee22f/IGVFFI8791HIOU.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8791HIOU.yaml"

echo "Decompressing IGVFFI1685FCPZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/bed8636d-e824-436e-b7f8-bdf6b4ca256c/IGVFFI1685FCPZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1685FCPZ.yaml"

echo "Decompressing IGVFFI0274PUWA.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/c08f8cf4-7988-4f5c-be2b-2e50c0c670a3/IGVFFI0274PUWA.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0274PUWA.yaml"

echo "Decompressing IGVFFI7933QCZK.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/c12c2df0-7028-4ecc-9ca9-ccf00f4cee02/IGVFFI7933QCZK.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7933QCZK.yaml"

echo "Decompressing IGVFFI5820TKFN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/ca02ec99-4754-407a-8b49-d64442b661c8/IGVFFI5820TKFN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5820TKFN.yaml"

echo "Decompressing IGVFFI6867QZGF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/cc8feeb9-597b-4b05-9020-98029ff1ee15/IGVFFI6867QZGF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6867QZGF.yaml"

echo "Decompressing IGVFFI0836FIUK.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/cf2d1128-7902-44c4-825c-66d7999f7fcf/IGVFFI0836FIUK.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0836FIUK.yaml"

echo "Decompressing IGVFFI4987DUPG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/e21f0698-1060-43c0-bcc9-c0607600e649/IGVFFI4987DUPG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4987DUPG.yaml"

echo "Decompressing IGVFFI1701AUAR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/e636cf6b-9458-466b-bf56-0609cdb814d2/IGVFFI1701AUAR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1701AUAR.yaml"

echo "Decompressing IGVFFI3074EJUE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/e934fdc2-8372-4748-8836-6c11936dcccb/IGVFFI3074EJUE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3074EJUE.yaml"

echo "Decompressing IGVFFI7947YZLX.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/ec232034-a21e-45f7-bb15-7a6e1a94fe0d/IGVFFI7947YZLX.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7947YZLX.yaml"

echo "Decompressing IGVFFI4446ALOB.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/f48fbe2c-c370-475a-90b9-bb816c28b156/IGVFFI4446ALOB.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4446ALOB.yaml"

echo "Decompressing IGVFFI2355JAYZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/f654cd20-16e3-4fcf-b4e8-fc01aa6b9585/IGVFFI2355JAYZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2355JAYZ.yaml"

echo "Decompressing IGVFFI4573QVNI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/fa9d96bf-30e8-4ad3-be8d-e1ecd34fd8b3/IGVFFI4573QVNI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4573QVNI.yaml"

echo "Decompressing IGVFFI9818AEJQ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_GEM-Xv3/2026_02_04/2025/08/06/fe0ba8f9-9e90-4cd0-8c5f-ce0f845698ea/IGVFFI9818AEJQ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9818AEJQ.yaml"

echo "Done\! Patched files in $PATCH_DIR"
gsutil ls "$PATCH_DIR/"
