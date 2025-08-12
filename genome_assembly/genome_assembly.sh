# commands used to generate genome assemblies from Oxford Nanopore Long read sequencing
# strain represents the BDSC strain numbers
# Author: Alex Samano, 2024


##### PREP READS #####
# convert basecalled reads to fastq
module load GCC/12.2.0 SAMtools/1.17

samtools bam2fq $strain"_ont_basecalled.bam" > $strain"_ont_basecalled.fastq"

# clean reads to remove Q<10 and length less than 10 kb
conda activate chopper_env
chopper -q 10 -l 10000 -i $strain"_ont_basecalled.fastq" > $strain"_ont_filt10kb.fastq"


##### ASSEMBLE #####

module load GCCcore/13.3.0 hifiasm/0.25.0

hifiasm --ont -o $strain"_hifiasm" -t 48 $strain"_ont_filt10kb.fastq"
# convert gfa to fasta
awk '/^S/{print ">"$2;print $3}' $strain"_hifiasm.bp.p_ctg.gfa" > $strain"_hifiasm_primary_contigs.fa"


##### DECONTAMINATE #####

# run BLAST
module load GCC/12.3.0  OpenMPI/4.1.5 BLAST+/2.14.1
export BLASTDB=/scratch/data/bio/blast/

blastn -query $strain"_hifiasm_primary_contigs.fa" -db nt -outfmt '6 qseqid staxids bitscore std' -num_threads 48 -max_target_seqs 1 -max_hsps 1  -evalue 1e-25 -out $strain.blast.out

# align back to contigs
minimap2 -ax map-ont  $strain"_hifiasm_primary_contigs.fa" $strain"_ont_filt10kb.fastq" | samtools sort -@48 -O BAM -o $strain"_reads2rawcontigs.bam" -
samtools index $strain"_reads2rawcontigs.bam"


# blob tools 
module load GCC/7.3.0-2.30  OpenMPI/3.1.1 BlobTools/20180528-Python-2.7.15
# make blob database
blobtools create -i $strain"_hifiasm_primary_contigs.fa" -b $strain"_reads2rawcontigs.bam" -t $strain.blast.out -o $strain
#make blob plots
blobtools view -i $strain.blobDB.json -o $strain"_blob"
blobtools plot -i $strain.blobDB.json 
# filter out microbial contigs
awk '$6 == "Arthropoda" || $6 == "no-hit"' $strain"_blob.$strain.blobDB.table.txt" | cut -f1 > $strain"_contigs2keep.txt"
blobtools seqfilter -i $strain"_hifiasm_primary_contigs.fa" -l $strain"_contigs2keep.txt"


##### POLISH #####
module load GCC/12.3.0  OpenMPI/4.1.5 medaka/1.12.0 
medaka_consensus -i $strain"_ont_filt10kb.fastq" -d $strain"_hifiasm_primary_contigs_clean.fa" -o medaka_wholeasm -t 24

##### SCAFFOLD #####

# mask repeats
module load GCC/11.3.0  OpenMPI/4.1.4 RepeatMasker/4.1.5
RepeatMasker -species 7227 -pa 24 $strain"_hifiasm_primary_contigs_clean_polished.fa"

module purge

# align masked contigs to masked reference (no Y because only females were sequenced)
module load GCC/11.2.0  OpenMPI/4.1.1  MUMmer/3.23
masked_ref=/scratch/user/asamano/dmel/ref/r6_49/dmel-r6.49-234X_masked.fa
nucmer -mum -prefix $strain"_masked_contigs2iso1" $masked_ref $strain"_hifiasm_primary_contigs_clean_polished.fa.masked"

# apply filters
delta-filter -1 $strain"_masked_contigs2iso1.delta" > $strain"_masked_contigs2iso1.1delta"
delta-filter -m $strain"_masked_contigs2iso1.delta" > $strain"_masked_contigs2iso1.mdelta"

# run mscaffolder (https://github.com/mahulchak/mscaffolder)
mscaffolder=/scratch/user/asamano/tools/mscaffolder/mscaffolder
$mscaffolder -d1 $strain"_masked_contigs2iso1.1delta" -md $strain"_masked_contigs2iso1.mdelta" -f $strain"_hifiasm_primary_contigs_clean_polished.fa" -ul n > $strain"_hifiasm_primary_scaffs.fa"

