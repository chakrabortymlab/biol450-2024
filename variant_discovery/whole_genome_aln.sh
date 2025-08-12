# commands to perform whole genome alignment and SV calling with SVMU

# Author: Alex Samano, 2025


##### NUCMER #####
contigs=$strain"_hifiasm_primary_contigs_clean_polished.fa"
scaffs=$strain"_hifiasm_primary_scaffs.fa"

#run nucmer
module load  GCCcore/13.2.0  MUMmer/4.0.0rc1
ref=/scratch/user/asamano/dmel/ref/r6_49/dmel_iso1_234X.fasta 

nucmer -p $strain"_scaffs_iso1ref" $ref $scaffs # align scaffolds
nucmer -p $strain"_contigs_iso1ref" $ref $contigs # align contigs


# get alignment coordinates
show-coords $strain"_scaffs_iso1ref.delta" > $strain"_scaffs_iso1ref_coords.txt"
show-coords $strain"_contigs_iso1ref.delta" > $strain"_contigs_iso1ref_coords.txt"


##### RUN SVMU #####
# available from https://github.com/mahulchak/svmu

svmu=/scratch/user/asamano/dmel/tools/svmu/svmu

#run svmu
$svmu $strain"_scaffs_iso1ref.delta" $ref $scaffs l null $strain"_scaffs_iso1ref" # scaffs
$svmu $strain"_contigs_iso1ref.delta" $ref $contigs l null $strain"_contigs_iso1ref" # contigs

# annotate repeats
module load GCC/11.3.0  OpenMPI/4.1.4 RepeatMasker/4.1.5

RepeatMasker -species 7227 -pa 24 -gff $scaffs
