# commands for calling variants from read alignment
# Author: Alex Samano, 2025

##### MAP READS #####
module load GCC/12.3.0 minimap2/2.26 SAMtools/1.18

ref=/scratch/user/asamano/dmel/ref/r6_49/dmel-all-chromosome-r6.49.fasta 
# map strain reads to iso1 reference then convert SAM file to sorted BAM file
minimap2 -t 48 -ax map-ont $ref $strain"_ont_basecalled.fastq" | samtools sort -@48 -O BAM -o $strain"_iso1ref_sorted.bam" -

# index BAM file
samtools index $strain"_iso1ref_sorted.bam"

##### CALL SMALL VARIANTS #####

# Pepper-Margin-DeepVariant Pipeline
# singularity image downloaded Jan 9, 2025
ref=/scratch/user/asamano/dmel/ref/r6_49/dmel-all-chromosome-r6.49.fasta 

singularity exec /sw/hprc/sw/containers/pepper/pepper_deepvariant_r0.8.sif \
run_pepper_margin_deepvariant call_variant \
-b $strain"_iso1ref_sorted.bam" \
-f $ref \
-o "pmdv_"$strain \
-p $strain"_pmdv" \
-t 24 \
--ont_r10_q20


##### CALL SVS #####
module load GCC/13.3.0  OpenMPI/5.0.3 Sniffles/2.6.1

sniffles -i $strain"_iso1ref_sorted.bam" -v $strain"_iso1ref_snifflesSV.vcf" --reference $ref

