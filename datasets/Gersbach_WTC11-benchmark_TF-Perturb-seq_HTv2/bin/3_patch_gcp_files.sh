#\!/bin/bash
set -e

# Patch files for Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2
PATCH_DIR="gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/patch"

echo "Decompressing IGVFFI9487JPEN.tsv"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2024/12/05/996cf65e-5b5c-494d-b8d3-f1f1e8b0181c/IGVFFI9487JPEN.tsv.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9487JPEN.tsv"

echo "Decompressing IGVFFI5765HMZH.tsv"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/04/cc3d00eb-118c-4230-945f-3cdb847cd849/IGVFFI5765HMZH.tsv.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5765HMZH.tsv"

echo "Decompressing IGVFFI0579IJOO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/03a1bc3f-2e06-4c4d-9566-20bbcc9e7f71/IGVFFI0579IJOO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0579IJOO.yaml"

echo "Decompressing IGVFFI0477KDVO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/0851c1fa-216d-4261-913a-989a682c4e7e/IGVFFI0477KDVO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0477KDVO.yaml"

echo "Decompressing IGVFFI4450UEXH.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/090bd2c0-330a-41a0-a741-76cc4d4ff06a/IGVFFI4450UEXH.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4450UEXH.yaml"

echo "Decompressing IGVFFI1373LVQN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/0c791803-bb31-45ca-9ec0-5963fa64953c/IGVFFI1373LVQN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1373LVQN.yaml"

echo "Decompressing IGVFFI6643CAVQ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/12be6d4b-a15f-464e-862e-47bcf67354a4/IGVFFI6643CAVQ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6643CAVQ.yaml"

echo "Decompressing IGVFFI0467ISAP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/12c38f42-adff-4714-91dd-baba98693861/IGVFFI0467ISAP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0467ISAP.yaml"

echo "Decompressing IGVFFI2257PVCP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/1485184e-bebb-4866-b12a-9ef57e250900/IGVFFI2257PVCP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2257PVCP.yaml"

echo "Decompressing IGVFFI5169TPHR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/1489b20a-3bde-4a8e-80e0-49cefa29b259/IGVFFI5169TPHR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5169TPHR.yaml"

echo "Decompressing IGVFFI5913GHNZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/17af4d9c-0f03-45dc-b4ee-ed1acbfb478d/IGVFFI5913GHNZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5913GHNZ.yaml"

echo "Decompressing IGVFFI4417TKFK.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/1ea04613-ba41-4f0c-b1bd-facb0a96bb0e/IGVFFI4417TKFK.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4417TKFK.yaml"

echo "Decompressing IGVFFI5680TKAM.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/2a1343c8-a78f-4a9a-81c9-de3052024492/IGVFFI5680TKAM.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5680TKAM.yaml"

echo "Decompressing IGVFFI5560PMPN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/2c62d63d-ab5d-4982-924b-cb187c5f178f/IGVFFI5560PMPN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5560PMPN.yaml"

echo "Decompressing IGVFFI3182NZQP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/350a03be-313c-4e21-a561-fbe063b4f63f/IGVFFI3182NZQP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3182NZQP.yaml"

echo "Decompressing IGVFFI9487ZLTF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/3fe74a7f-cd6e-4d9a-949d-2b225c521bab/IGVFFI9487ZLTF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9487ZLTF.yaml"

echo "Decompressing IGVFFI8970YTOV.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/418cca94-cd20-4587-80c3-a22c4d411d1c/IGVFFI8970YTOV.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8970YTOV.yaml"

echo "Decompressing IGVFFI0199SLRN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/4871d7ee-0cb3-4f9e-b27c-1ff947a36151/IGVFFI0199SLRN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0199SLRN.yaml"

echo "Decompressing IGVFFI0936QREW.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/4e234348-d6ae-414a-88c2-f378b4390f55/IGVFFI0936QREW.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0936QREW.yaml"

echo "Decompressing IGVFFI2534RGZE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/4e4cb2f8-3af7-4b3f-88eb-617f0129293d/IGVFFI2534RGZE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2534RGZE.yaml"

echo "Decompressing IGVFFI9746QYCG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/515f9fdc-facf-461f-bead-94345f25694b/IGVFFI9746QYCG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9746QYCG.yaml"

echo "Decompressing IGVFFI9593AAFI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/5164e95d-491c-4dc6-8243-ecc906b357e6/IGVFFI9593AAFI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9593AAFI.yaml"

echo "Decompressing IGVFFI6648OTPT.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/53a38fa9-d626-48e8-a111-034fc4cc8882/IGVFFI6648OTPT.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI6648OTPT.yaml"

echo "Decompressing IGVFFI5900ANFB.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/5914d5f3-b0f4-4b8f-84c1-3e43aad6468d/IGVFFI5900ANFB.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5900ANFB.yaml"

echo "Decompressing IGVFFI9684UINN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/5ec6a4ea-90ad-4171-b7fa-74bd8e3f9a3c/IGVFFI9684UINN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9684UINN.yaml"

echo "Decompressing IGVFFI8291UJNI.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/60320ae7-4a98-4091-b717-c42e37219af3/IGVFFI8291UJNI.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8291UJNI.yaml"

echo "Decompressing IGVFFI2755QRLF.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/6498ea83-8145-4215-b62d-6f60a6345d90/IGVFFI2755QRLF.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2755QRLF.yaml"

echo "Decompressing IGVFFI8573ECKZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/65ed695b-7492-4757-9ef7-97673fafdce4/IGVFFI8573ECKZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8573ECKZ.yaml"

echo "Decompressing IGVFFI9127TSTV.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/661ecdbb-59ea-4cd2-b169-245870f447ae/IGVFFI9127TSTV.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9127TSTV.yaml"

echo "Decompressing IGVFFI1957ERIA.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/6cbfca7a-1de7-4ec7-b835-60a1a6b31483/IGVFFI1957ERIA.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1957ERIA.yaml"

echo "Decompressing IGVFFI9735LVDP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/71483adb-9c52-43d9-a3e1-5302fb557660/IGVFFI9735LVDP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9735LVDP.yaml"

echo "Decompressing IGVFFI4556HOVD.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/7a2c6868-0981-4080-8946-a29faa538629/IGVFFI4556HOVD.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4556HOVD.yaml"

echo "Decompressing IGVFFI9143GPKP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/7d1a8ce4-6063-465a-afd9-d4674b2416da/IGVFFI9143GPKP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9143GPKP.yaml"

echo "Decompressing IGVFFI9979FCKE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/7e7bba44-2443-4a5c-bb9b-b7bbff113ccd/IGVFFI9979FCKE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9979FCKE.yaml"

echo "Decompressing IGVFFI4189ERPQ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/7fd1566a-448d-4a7c-b813-a08898faa525/IGVFFI4189ERPQ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4189ERPQ.yaml"

echo "Decompressing IGVFFI1490SRKG.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/88b18826-1330-45ff-82ed-6bc054b729af/IGVFFI1490SRKG.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1490SRKG.yaml"

echo "Decompressing IGVFFI0872WZAZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/8c2e6c6b-18fe-4671-954e-7b78b8c4d4b1/IGVFFI0872WZAZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0872WZAZ.yaml"

echo "Decompressing IGVFFI0539DHAP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/8f0f5eda-04e8-46c3-a175-880bb641cefa/IGVFFI0539DHAP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0539DHAP.yaml"

echo "Decompressing IGVFFI7632IZCM.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/90221420-9a22-4a7d-9107-4ee6d6aff6b3/IGVFFI7632IZCM.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7632IZCM.yaml"

echo "Decompressing IGVFFI9626QPYJ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/946a0c69-c9b0-4244-b78f-78428de2789f/IGVFFI9626QPYJ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9626QPYJ.yaml"

echo "Decompressing IGVFFI5476ZWRC.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/94e9a07e-6e9f-4970-ad61-431fc2b41410/IGVFFI5476ZWRC.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5476ZWRC.yaml"

echo "Decompressing IGVFFI3147SLSZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/9a67e1c4-a334-486e-98d1-5fb757887f8e/IGVFFI3147SLSZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI3147SLSZ.yaml"

echo "Decompressing IGVFFI2756QIPC.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/9b64c5ac-53f2-4dfe-8bd7-00537e2d5169/IGVFFI2756QIPC.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2756QIPC.yaml"

echo "Decompressing IGVFFI7485FQKE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/9f8ff998-1859-4045-926f-a39e83889c2c/IGVFFI7485FQKE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7485FQKE.yaml"

echo "Decompressing IGVFFI9789XKYX.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/a17b9d54-f902-41f0-93a2-b73cfe381b35/IGVFFI9789XKYX.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9789XKYX.yaml"

echo "Decompressing IGVFFI1690RTSC.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/a1c56e1d-3f96-4515-a626-242346acaa3d/IGVFFI1690RTSC.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1690RTSC.yaml"

echo "Decompressing IGVFFI2036NNRP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/a7698d62-f4ae-4f7f-8f49-a327580b0fcd/IGVFFI2036NNRP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2036NNRP.yaml"

echo "Decompressing IGVFFI2767XXUK.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/a947aba7-49f2-4fa1-8f92-77f2f8b7a048/IGVFFI2767XXUK.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2767XXUK.yaml"

echo "Decompressing IGVFFI9467YPMD.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/aa7c7852-b07e-4b05-81d7-79eb6657a92c/IGVFFI9467YPMD.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9467YPMD.yaml"

echo "Decompressing IGVFFI8944YQAH.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/af3af416-0f59-4484-ae69-b763ca8ec3f6/IGVFFI8944YQAH.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI8944YQAH.yaml"

echo "Decompressing IGVFFI7807YOBM.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/b3aa027d-fd5d-4399-9b82-5def616bfebf/IGVFFI7807YOBM.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7807YOBM.yaml"

echo "Decompressing IGVFFI2294CBIY.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/b8242412-54b6-43a0-afda-d1b75fe8a33b/IGVFFI2294CBIY.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2294CBIY.yaml"

echo "Decompressing IGVFFI4368QZBJ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/ce776c08-e9cd-4ffe-a168-65cf2521c488/IGVFFI4368QZBJ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI4368QZBJ.yaml"

echo "Decompressing IGVFFI5578YVQU.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/d6740acf-0e08-405d-a1d1-b292d3a22186/IGVFFI5578YVQU.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5578YVQU.yaml"

echo "Decompressing IGVFFI9946JRCZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/d6a003a6-6e56-421d-8aa4-d294fb8a65e3/IGVFFI9946JRCZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9946JRCZ.yaml"

echo "Decompressing IGVFFI5716EDJN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/d920d638-f6c0-47df-9c39-ecfe6a045218/IGVFFI5716EDJN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5716EDJN.yaml"

echo "Decompressing IGVFFI9006WHIZ.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/dd46539e-4f21-4652-9cc4-7f7fbb93b9df/IGVFFI9006WHIZ.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9006WHIZ.yaml"

echo "Decompressing IGVFFI7426GTQS.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/ddb71812-5da7-472d-96c1-6120d03ccb15/IGVFFI7426GTQS.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7426GTQS.yaml"

echo "Decompressing IGVFFI9501ZDCE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/e041c609-f1f3-4efc-8601-349648538379/IGVFFI9501ZDCE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9501ZDCE.yaml"

echo "Decompressing IGVFFI7643CUNP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/e04f5652-f56b-47eb-b165-6e7062ff81fa/IGVFFI7643CUNP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7643CUNP.yaml"

echo "Decompressing IGVFFI0326NMJT.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/e19567bd-c444-4b37-8ae5-0a6b6d012a2c/IGVFFI0326NMJT.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI0326NMJT.yaml"

echo "Decompressing IGVFFI9254GKSO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/e5702a97-1c36-4116-9c55-1bdc76317e73/IGVFFI9254GKSO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9254GKSO.yaml"

echo "Decompressing IGVFFI2148JIGE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/e9bb6f6c-0b5c-4c97-af49-ded5790a8165/IGVFFI2148JIGE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2148JIGE.yaml"

echo "Decompressing IGVFFI2009DPQO.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/ed4b7945-83aa-4a74-8b3f-52ce83d33dc3/IGVFFI2009DPQO.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2009DPQO.yaml"

echo "Decompressing IGVFFI2576JLSC.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/edae55eb-ddc4-43a1-bae7-c745e17737a0/IGVFFI2576JLSC.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI2576JLSC.yaml"

echo "Decompressing IGVFFI9026VIVR.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/f0b1918b-9cc4-48d8-80fc-113d10f71379/IGVFFI9026VIVR.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI9026VIVR.yaml"

echo "Decompressing IGVFFI1572MKFP.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/f2d95ec6-0919-49f9-b5d7-c7b13a29082e/IGVFFI1572MKFP.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI1572MKFP.yaml"

echo "Decompressing IGVFFI5953SIMN.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/f67f18f4-2e80-44ad-9a94-a7d2ea14157b/IGVFFI5953SIMN.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5953SIMN.yaml"

echo "Decompressing IGVFFI7449MVEE.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/f9aa3fe1-955c-42c9-9602-0b05e4de11a9/IGVFFI7449MVEE.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI7449MVEE.yaml"

echo "Decompressing IGVFFI5554CQTU.yaml"
gsutil cat "gs://igvf-pertub-seq-pipeline-data/Gersbach_WTC11-benchmark_TF-Perturb-seq_HTv2/2026_02_04/2025/08/06/fc508ab6-7557-4995-863f-3ede921eff77/IGVFFI5554CQTU.yaml.gz" | gunzip | gsutil cp - "$PATCH_DIR/IGVFFI5554CQTU.yaml"

echo "Done\! Patched files in $PATCH_DIR"
gsutil ls "$PATCH_DIR/"
