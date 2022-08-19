# Generating Reciprocal Best Structural Hits

## Requirements

Foldseek (https://github.com/steineggerlab/foldseek) to be installed in /bin of the project root <br>
Python >= 3.9 <br>
Config file with completed fields

## Usage
`python find_RBSH.py config.ini --runfoldseek [yes, no]`

# Generating Reciprocal Best Hits

## Requirements
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (version  2.12.0+) <br>
Python >= 3.9

## Usage
`python3.9 find_RBH.py --fasta_A proteome_seqs_A.fasta --fasta_B proteome_seqs_B.fasta --output_file RBH_A_B.csv`