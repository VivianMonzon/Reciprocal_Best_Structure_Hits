#!/bin/bash                                                                                                                                                                                                 

declare -A proteome_ids
proteome_ids["H_sapiens"]="H_sapiens_UP000005640_9606_HUMAN_v2"
proteome_ids["D_melanogaster"]="D_melanogaster_UP000000803_7227_DROME_v2"
proteome_ids["C_elegans"]="C_elegans_UP000001940_6239_CAEEL_v2"
proteome_ids["S_cerevisiae"]="S_cerevisiae_UP000002311_559292_YEAST_v2"
proteome_ids["S_pombe"]="S_pombe_UP000002485_284812_SCHPO_v2"

main(){
    preparations
    remove_long_proteins
    get_AF_structure_ids
    unzip_structure_files
    get_seq_from_cif_file
    domain_analysis
    iupred_analysis
    rm_low_quality_regions
}

preparations(){
    mkdir data
    mkdir structure_db
    for var in 'C_elegans' 'D_melanogaster' 'S_cerevisiae' 'S_pombe' 'H_sapiens'
    do
	mkdir structure_db/"${proteome_ids["$var"]}"/
    done
    wget -O structure_db/D_melanogaster_UP000000803_7227_DROME_v2.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000803_7227_DROME_v3.tar
    wget -O structure_db/C_elegans_UP000001940_6239_CAEEL_v2.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000001940_6239_CAEEL_v2.tar
    wget -O structure_db/S_cerevisiae_UP000002311_559292_YEAST_v2.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v2.tar
    wget -O structure_db/S_pombe_UP000002485_284812_SCHPO_v2.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002485_284812_SCHPO_v2.tar
    wget -O structure_db/H_sapiens_UP000005640_9606_HUMAN_v2.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v2.tar

    for var in 'C_elegans' 'D_melanogaster' 'S_cerevisiae' 'S_pombe' 'H_sapiens'
    do
        tar -xf structure_db/"${proteome_ids["$var"]}".tar -C structure_db/"${proteome_ids["$var"]}"/
    done
}

remove_long_proteins(){
    mkdir structure_db/H_sapiens_wo_long_proteins/
    ls structure_db/H_sapiens_UP000005640_9606_HUMAN_v2/*.pdb.gz > structure_db/H_sapiens_AF_ids_w_long_proteins.csv
    python3.9 remove_long_proteins.py
}

get_AF_structure_ids(){
    for var in 'C_elegans' 'D_melanogaster' 'S_cerevisiae' 'S_pombe'
    do
	ls structure_db/"${proteome_ids["$var"]}"/*.pdb.gz | cut -d'-' -f2 > structure_db/${var}_AF_ids.csv
    done
    ls structure_db/H_sapiens_wo_long_proteins/*.pdb.gz | cut -d'-' -f2 > structure_db/H_sapiens_AF_ids_wo_long_prots.csv
}

unzip_structure_files(){
    for i in structure_db/S_cerevisiae_UP000002311_559292_YEAST_v2/*.gz; do gunzip $i; done
    for i in structure_db/C_elegans_UP000001940_6239_CAEEL_v2/*.gz; do gunzip $i; done
    for i in structure_db/D_melanogaster_UP000000803_7227_DROME_v2/*.gz; do gunzip $i; done
    for i in structure_db/H_sapiens_wo_long_proteins/*.gz; do gunzip $i; done
    for i in structure_db/S_pombe_UP000002485_284812_SCHPO_v2/*.gz; do gunzip $i; done
}

get_seq_from_cif_file(){
    for i in structure_db/S_cerevisiae_UP000002311_559292_YEAST_v2/*.cif; do \
	python3.9 extrac_seqs.py --cif_file $i --out_seq_fa data/S_cerevisiae_AF_cif_seqs.fasta; done
    for i in structure_db/C_elegans_UP000001940_6239_CAEEL_v2/*.cif; do \
	python3.9 extrac_seqs.py --cif_file $i --out_seq_fa data/C_elegans_AF_cif_seqs.fasta; done
    for i in structure_db/D_melanogaster_UP000000803_7227_DROME_v2/*.cif; do \
	python3.9 extrac_seqs.py --cif_file $i --out_seq_fa data/D_melanogaster_AF_cif_seqs.fasta; done
    for i in structure_db/H_sapiens_wo_long_proteins/*.cif; do \
	python3.9 extrac_seqs.py --cif_file $i --out_seq_fa data/H_sapiens_AF_cif_seqs.fasta; done
    for i in structure_db/S_pombe_UP000002485_284812_SCHPO_v2/*.cif; do \
	python3.9 extrac_seqs.py --cif_file $i --out_seq_fa data/S_pombe_AF_cif_seqs.fasta; done
}

domain_analysis(){
    for var in 'C_elegans' 'D_melanogaster' 'H_sapiens' 'S_cerevisiae' 'S_pombe'
    do
        hmmsearch --cut_ga --domtblout data/Pfam_doms_${var}.tbl Pfam-A_35.hmm data/${var}_AF_cif_seqs.fasta > data/Pfam_doms_${var}.out
    done
}

iupred_analysis(){
    mkdir data/iupred_results/
    for var in 'C_elegans' 'D_melanogaster' 'S_cerevisiae' 'S_pombe' 'H_sapiens'
    do
        python3.9 iupred2a.py data/${var}_AF_cif_seqs.fasta 'long' > data/iupred_results/${var}_AF_cif_iupred.csv
	python3.9 iupred_adapt.py --fh_in data/iupred_results/${var}_AF_cif_iupred.csv \
		  --fh_out data/iupred_results/${var}_AF_cif_iupred_adapt.csv
    done
}

rm_low_quality_regions(){
    # S. pombe:
    mkdir S_pombe_structures_adapted/
    for i in structure_db/S_pombe_UP000002485_284812_SCHPO_v2/*.pdb; do \
        python3.9 rm_low_plddt_regions.py --pdb_in $i \
                  --pdb_out S_pombe_structures_adapted/$(basename $i .pdb)_cut.pdb; done
    # S. cerevisiae:
    mkdir S_cerevisiae_structures_adapted/
    for i in structure_db/S_cerevisiae_UP000002311_559292_YEAST_v2/*.pdb; do \
        python3.9 rm_low_plddt_regions.py --pdb_in $i \
                  --pdb_out S_cerevisiae_structures_adapted/$(basename $i .pdb)_cut.pdb; done
    # C. elegans:
    mkdir C_elegans_structures_adapted/
    for i in structure_db/C_elegans_UP000001940_6239_CAEEL_v2/*.pdb; do \
        python3.9 rm_low_plddt_regions.py --pdb_in $i \
                  --pdb_out C_elegans_structures_adapted/$(basename $i .pdb)_cut.pdb; done
    # H. spaiens:
    mkdir H_sapiens_structures_adapted/
    for i in structure_db/H_sapiens_wo_long_proteins/*.pdb; do \
        python3.9 rm_low_plddt_regions.py --pdb_in $i \
                  --pdb_out H_sapiens_structures_adapted/$(basename $i .pdb)_cut.pdb; done
    # D. melanogaster:
    mkdir D_melanogaster_structures_adapted/
    for i in structure_db/D_melanogaster_AF/*.pdb; do \
        python3.9 rm_low_plddt_regions.py --pdb_in $i \
                  --pdb_out D_melanogaster_structures_adapted/$(basename $i .pdb)_cut.pdb; done
}

main
