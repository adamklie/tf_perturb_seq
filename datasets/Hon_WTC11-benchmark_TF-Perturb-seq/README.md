# Hon_WTC11-benchmark_TF-Perturb-seq

## TODO
- [ ] Update guide metadata to latest spec
- [ ] Rerun through pipeline 
- [ ] Get CellRanger and PySpade output from Mpathi
- [ ]

## Experimental design

| **System** | **Sequencing Tool** | **Lab** | **Contact** | **Technology (3'  , 5')** | **MOI  (HIGH or LOW)** | **Number of cells sequenced** |
| --- | --- | --- | --- | --- | --- | --- |
| iPSC (WTC11) | 10x 5' (w/HTO) | Hon | Mpathi | 5' to 3' | HIGH | 750K |

## Data portal 
[IGVFDS4761PYUO](https://data.igvf.org/analysis-sets/IGVFDS4761PYUO/)

## Nextflow samplesheet
This tooka few iterations to get right. I used Ian's initial download scripts to download the raw fastqs and relevant metadata and upload to GCP (took some iteration). Seqpsecs were not on the portal yet so I grabbed the generic ones from the pipeline repo and uploaded them to GCP manually. There were a couple other manual updates I had to do as there was some funkiness with the first row of the samplesheet. The version actually used can be found here: https://www.synapse.org/Synapse:syn70753575

## Guide metadata
Generated v1 metadata as described [here](https://docs.google.com/presentation/d/1uqnACxqehcuuI2o32ydW8EsVTQ-3j6oiuTxBTPwJM8I/edit?slide=id.p#slide=id.p)
Recent edits to metadata spec require that this be updated and the pipeline rerun.

## Uniform processing and inference
Generated v1 pipeline outputs on GCP. I think the html report failed to generate so the relevant outputs can only be found in the work directory. This is only really something you can navigate if you have access to all the logs of the nextflow run via Seqera Tower (I believe Lucas has this access).

Stored the inference MuData copied over here: https://www.synapse.org/Synapse:syn70753405

QC Analysis slides: https://docs.google.com/presentation/d/1E5aUUtMdyWWutOi9Uj24nlXjYkLPOfgm7N2XbEqjrtE/edit?slide=id.g39e356c57dd_0_11#slide=id.g39e356c57dd_0_11

## Gene program discovery
https://www.synapse.org/Synapse:syn70776514
