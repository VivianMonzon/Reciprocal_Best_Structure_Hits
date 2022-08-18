import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fh_foldseek')
parser.add_argument('--fh_reci_blastp')
parser.add_argument('--blastp_all_p1')
parser.add_argument('--blastp_all_p2')
parser.add_argument('--all_foldseek_results')
parser.add_argument('--common')
parser.add_argument('--only_foldseek_out')
parser.add_argument('--only_blastp_foldseek_all')
parser.add_argument('--only_blastp_foldseek_cov_thres')
args = parser.parse_args()

cols = ['prot1', 'prot2', 'qlen_p1', 'tlen_p1', 'fident_p1', 'alnlen_p1',
        'mismatch_p1', 'qstart_p1', 'qend_p1', 'tstart_p1', 'tend_p1',
        'evalue_p1', 'bits_p1', 'qlen_p2', 'tlen_p2', 'fident_p2', 'alnlen_p2',
        'mismatch_p2', 'qstart_p2', 'qend_p2', 'tstart_p2', 'tend_p2',
        'evalue_p2', 'bits_p2']
df_foldseek = pd.read_csv(args.fh_foldseek, names=cols)

df_foldseek['qlen_p1'] = df_foldseek['qlen_p1'].apply(
    lambda x: int(x.split("['")[1]))
df_foldseek['bits_p1'] = df_foldseek['bits_p1'].apply(
    lambda x: int(x.split("']")[0]))
df_foldseek['qlen_p2'] = df_foldseek['qlen_p2'].apply(
    lambda x: int(x.split("['")[1]))
df_foldseek['bits_p2'] = df_foldseek['bits_p2'].apply(
    lambda x: int(x.split("']")[0]))

df_foldseek['q_cov_p1'] = (
    (df_foldseek['qend_p1'] - df_foldseek['qstart_p1'])/df_foldseek['qlen_p1'])*100
df_foldseek['q_cov_p2'] = (
    (df_foldseek['qend_p2'] - df_foldseek['qstart_p2'])/df_foldseek['qlen_p2'])*100
df_foldseek.to_csv(args.all_foldseek_results, index=False)

df_blastp = pd.read_csv(args.fh_reci_blastp)
df_blastp_cols = df_blastp.columns
df_blastp['reciprocals'] = df_blastp[df_blastp.columns[0]] + '_' + \
    df_blastp[df_blastp.columns[1]]
blastp_recis = df_blastp['reciprocals'].tolist()
print('Number RBH hits:')
print(len(blastp_recis))

df_foldseek_all = df_foldseek
df_foldseek_all['reciprocals'] = df_foldseek_all['prot1'] + '_' + \
    df_foldseek_all['prot2']
foldseek_recis_all = df_foldseek_all['reciprocals'].tolist()
print('Number RBSH all hits (no cov. threshold):')
print(len(foldseek_recis_all))

same_recis_all = [x for x in foldseek_recis_all if x in blastp_recis]
print('Common hits (foldseek all - no cov. threshold):')
print(len(same_recis_all))

only_foldseek_all = [x for x in foldseek_recis_all if x not in blastp_recis]
print('Only RBSH all (no RBSH cov. threshold):')
print(len(only_foldseek_all))

only_blastp_no_foldseek_all = [
    x for x in blastp_recis if x not in foldseek_recis_all]
print('Only RBH (no RBSH cov. thres.):')
print(len(only_blastp_no_foldseek_all))

fh_blastp_out = open(args.only_blastp_foldseek_all, 'w')
for x in only_blastp_no_foldseek_all:
    fh_blastp_out.write('{}\n'.format(x))

print('RBSH with 75% coverage threshold:')
df_foldseek = df_foldseek[(df_foldseek['q_cov_p1'] > 75 ) & (
    df_foldseek['q_cov_p2'] > 75)]
df_foldseek['reciprocals'] = df_foldseek['prot1'] + '_' + df_foldseek['prot2']

foldseek_recis = df_foldseek['reciprocals'].tolist()
print('RBSH hits:')
print(len(foldseek_recis))

same_recis = [x for x in foldseek_recis if x in blastp_recis]
print('Common hits:')
print(len(same_recis))
common_fh = open(args.common, 'w')
for x in same_recis:
    common_fh.write('{}\n'.format(x))
 
only_foldseek = [x for x in foldseek_recis if x not in blastp_recis]
print('Only RBSH:')
print(len(only_foldseek))

only_blastp = [x for x in blastp_recis if x not in foldseek_recis]
print('Only RBH:')
print(len(only_blastp))

fh_blastp_out_fs_cov = open(args.only_blastp_foldseek_cov_thres, 'w')
for x in only_blastp:
    fh_blastp_out_fs_cov.write('{}\n'.format(x))

blastp_Fields = ['query', 'subject', 'identity', 'alignment length',
                 'mismatches', 'gap opens', 'q. start', 'q. end',
                 'query length', 's. start', 's. end', 'subject length',
                 'evalue', 'bit score', 'subject coverage']

all_blastp_prot1 = pd.read_csv(args.blastp_all_p1, sep='\t',
                               names=blastp_Fields)
all_blastp_prot1['p1_query_p1_p2'] = all_blastp_prot1['query'] + '_' + \
    all_blastp_prot1['subject']
p1_query_p1_p2 = list(set(all_blastp_prot1['p1_query_p1_p2'].tolist()))
all_blastp_prot2 = pd.read_csv(args.blastp_all_p2, sep='\t',
                               names=blastp_Fields)
all_blastp_prot2['p2_query_p2_p1'] = all_blastp_prot2['subject'] + '_' + \
    all_blastp_prot2['query']
p2_query_p2_p1 = list(set(all_blastp_prot2['p2_query_p2_p1'].tolist()))
all_blastp_reciprocal = [x for x in p1_query_p1_p2 if x in p2_query_p2_p1]

print('Only RBSH - but found with RBH with a lower score or coverage:')
print(len([x for x in only_foldseek if x in all_blastp_reciprocal]))

print('Only RBSH - not found with RBH with a lower score or coverage:')
only_fodlseek_no_blastp_at_all = [
    x for x in only_foldseek if x not in all_blastp_reciprocal]
print(len(only_fodlseek_no_blastp_at_all))

fh_out = open(args.only_foldseek_out, 'w')
for x in only_fodlseek_no_blastp_at_all:
    fh_out.write('{}\n'.format(x))
