# commands for quality assessment of genomes
# Author: Alex Samano, 2024

##### QUAST #####

module load GCC/9.3.0  OpenMPI/4.0.3 QUAST/5.0.2-Python-3.8.2

quast -o $strain"_contigs_quast" $strain"_hifiasm_primary_contigs_clean_polished.fa"
quast -o $strain"_scaffs_quast" $strain"_hifiasm_primary_scaffs.fa"

##### BUSCO #####
module load GCC/12.2.0  OpenMPI/4.1.4 BUSCO/5.7.1

# diptera database downloaded July 8,2024
diptera_db=/scratch/group/chakraborty_lab/alex/diptera_odb10
busco -i $strain"_hifiasm_primary_contigs_clean_polished.fa" -m genome -l $diptera_db -c 24 -o $strain"_dip_busco" --offline

