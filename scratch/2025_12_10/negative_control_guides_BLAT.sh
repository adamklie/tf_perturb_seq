#mamba activate pipelines

genome=/cellar/users/aklie/data/ref/genomes/hg38/hg38.fa
guides=/cellar/users/aklie/data/datasets/tf_perturb_seq/scratch/2025_12_10/negative_control_guides.fasta
output=/cellar/users/aklie/data/datasets/tf_perturb_seq/scratch/2025_12_10/2025_12_10_negative_control_BLAT_v2.psl

cmd="blat -t=dna -q=dna $genome $guides -tileSize=10 -stepSize=1 -minScore=19 -minMatch=1 -minIdentity=95 -noSimpRepMask $output"
echo $cmd
eval $cmd
