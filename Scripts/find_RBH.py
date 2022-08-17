import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_A', required=True,
                    help='Sequences of organism A in Fasta format')
parser.add_argument('--fasta_B', required=True,
                    help='Sequences of organism B in Fasta format')
parser.add_argument('--output_file', required=True,
                    help='csv output of RBH hits')
args = parser.parse_args()

Fields = ['query', 'subject', 'identity', 'alignment length', 'mismatches',
          'gap opens', 'q. start', 'q. end', 'query length', 's. start',
          's. end', 'subject length', 'evalue', 'bit score',
          'subject coverage']


def run_blastp(seqs_query, seqs_db, out_fh, number):
    print('BLASTP search {}:'.format(number))
    os.system('makeblastdb -in {} -dbtype prot'.format(seqs_db))
    os.system(('blastp -db {} -query {} -out {} '
               '-outfmt "7 qseqid sseqid pident length mismatch gapopen '
               'qstart qend qlen sstart send slen evalue bitscore qcovs" '
               '-evalue 0.01'.format(seqs_db, seqs_query, out_fh)))

    
def find_best_hits(fh_blast_out):
    df = pd.read_csv(fh_blast_out, sep='\t', names=Fields, comment='#')
    df = df[df['evalue'].astype(float) <= 0.01]
    df = df[df['subject coverage'].astype(float) > 75]
    query_ids = list(set(df['query'].tolist()))
    df_top_hits = pd.DataFrame(columns=Fields, dtype=object)
    for q in query_ids:
        df_q = df[df['query'] == q]
        df_q_max = df_q[df_q['evalue'] == df_q['evalue'].min()]
        if df_q_max.shape[0] > 1:
            df_q_max = df_q_max[df_q_max['bit score'] == df_q_max['bit score'].max()]
        else:
            df_q_max = df_q_max
        df_top_hits = df_top_hits.append(df_q_max)
    return df_top_hits


run_blastp(args.fasta_A, args.fasta_B, 'blastp_B_db_vs_A_query.tsv',
           'A against B')
run_blastp(args.fasta_B, args.fasta_A, 'blastp_A_db_vs_B_query.tsv',
           'B against A')

df_one = find_best_hits('blastp_B_db_vs_A_query.tsv')
df_one = df_one.rename({'query': 'organism_A', 'subject': 'organism_B'},
                       axis=1)
df_one = df_one[['organism_A', 'organism_B']]

df_two = find_best_hits('blastp_A_db_vs_B_query.tsv')
df_two = df_two.rename({'query': 'organism_B', 'subject': 'organism_A'},
                       axis=1)
df_two = df_two[['organism_B', 'organism_A']]

df_RBH = df_one.merge(df_two, on=['organism_A', 'organism_B'],
                      how='inner')
df_RBH = df_RBH.drop_duplicates()
df_RBH.to_csv(args.output_file, index=False)
