# commands to create a pangenome graph from 11 scaffolded assemblies
# Author: Alex Samano, 2025

##### RUN MINIGRAPH-CACTUS #####

# run mgcactus on scaffolded assemblies and ISO1 chroms 2,3,X
# singularity image downloaded April 11, 2025
singularity exec /sw/hprc/sw/bio/containers/cactus/cactus_v2.5.2-gpu.sif cactus-pangenome job_storeClip biol450_genomes.txt \
--outDir clip_out_biol450_genomes_234X \ 
--outName biol450_genomes \
--reference ISO1 \
--vcf clip \
--gfa clip \
--gbz clip

# split into biallelic 
# shift to represent snps
bcftools norm -m- -f dmel-r6.49-234X.fa -o biol450_genomes_biallelic_norm.vcf biol450_genomes.vcf

#keep only indels
vcftools --keep-only-indels --recode --recode-INFO-all --out biol450_genomes_norm_INDELS.vcf --vcf biol450_genomes_biallelic_norm.vcf
