import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('--panther_p1', required=True)
parser.add_argument('--panther_p2', required=True)
parser.add_argument('--all_foldseek_results', required=True)
parser.add_argument('--blastp_all_p1', required=True)
parser.add_argument('--blastp_all_p2', required=True)
parser.add_argument('--png_out', required=True)
args = parser.parse_args()


def get_plot_values(FP_TP_list):
    y_values = [0]
    x_values = [0]
    y = 0
    x = 0
    for P in FP_TP_list:
        if P == 'TP':
            y += 1
            x += 0
            y_values.append(y)
            x_values.append(x)
        elif P == 'FP':
            y += 0
            x += 1
            y_values.append(y)
            x_values.append(x)
    y_values_series = pd.Series(y_values)
    x_values_series = pd.Series(x_values)
    return y_values_series, x_values_series
    

def get_TP_FP(df):
    df = df.merge(df_panther_p1, on='prot1', how='left').merge(df_panther_p2, on='prot2', how='left')
    df['Type'] = np.where(df['Family_p1'] == df['Family_p2'], 'TP', 'FP')
    TP = df[df['Type'] == 'TP'].shape[0]
    FP = df[df['Type'] == 'FP'].shape[0]
    Precision = round((TP/(TP+FP)), 3)
    print('Precision: {}'.format(Precision))
    df['mean_evalue'] = (df['evalue_p1'] + df['evalue_p2']) / 2
    df = df.sort_values(by='mean_evalue', ascending=True)
    y_values, x_values = get_plot_values(df['Type'].tolist())
    return y_values, x_values
        

def adapt_foldseek_results(fh_in):
    df = pd.read_csv(fh_in)
    df = df[['prot1', 'prot2', 'evalue_p1', 'evalue_p2',
             'q_cov_p1', 'q_cov_p2']]
    df['pairs'] = df['prot1'] + '_' + df['prot2']
    df_cov = df[(df['q_cov_p1'] > 75 ) & (df['q_cov_p2'] > 75)]
    return df, df_cov


def read_panther(fh_in, protein):
    df = pd.read_csv(fh_in, names=['ID', 'Subfamily', 'description'], sep='\t')
    df['Family'] = df['Subfamily'].apply(lambda x: x.split(':')[0])
    df = df.rename({'ID': 'prot{}'.format(protein),
                    'Family': 'Family_p{}'.format(protein),
                    'Subfamily': 'Subfamily_p{}'.format(protein),
                    'description': 'description_p{}'.format(protein)}, axis=1)
    df = df[['prot{}'.format(protein), 'Family_p{}'.format(protein)]]
    return df


Fields = ['query', 'subject', 'identity', 'alignment length', 'mismatches',
          'gap opens', 'q. start', 'q. end', 'query length', 's. start',
          's. end', 'subject length', 'evalue', 'bit score',
          'subject coverage']


def best_hit_blastp(fh_in, query_name, subject_name):
    df = pd.read_csv(fh_in, sep='\t', names=Fields, comment='#')
    df = df[df['evalue'] <= 0.01]
    df = df[df['subject coverage'] > 75]
    query_ids = list(set(df['query'].tolist()))
    df_top_hits = pd.DataFrame(columns=Fields, dtype=object)
    for q in query_ids:
        df_q = df[df['query'] == q]
        df_q_max = df_q[df_q['evalue'] == df_q['evalue'].min()]
        df_top_hits = df_top_hits.append(df_q_max)
    df_top_hits = df_top_hits[['query', 'subject', 'evalue']]
    df_top_hits = df_top_hits.drop_duplicates()
    df_top_hits = df_top_hits.rename({'query': query_name,
                                      'subject': subject_name}, axis=1)
    return df_top_hits


def read_blastp(fh_in_p1, fh_in_p2):
    df_p1 = best_hit_blastp(fh_in_p1, 'prot1', 'prot2')
    df_p2 = best_hit_blastp(fh_in_p2, 'prot2', 'prot1')
    df_reciprocal = df_p1.merge(df_p2, on=['prot1', 'prot2'], how='inner')
    df_reciprocal = df_reciprocal.rename(
        {'evalue_x': 'evalue_p1', 'evalue_y': 'evalue_p2'}, axis=1)
    return df_reciprocal


df_foldseek, df_foldseek_cov = adapt_foldseek_results(
    args.all_foldseek_results)

df_blastp = read_blastp(args.blastp_all_p1, args.blastp_all_p2)

df_panther_p1 = read_panther(args.panther_p1, 1)
df_panther_p2 = read_panther(args.panther_p2, 2)

print('Foldseek all:')
y_foldseek, x_foldseek = get_TP_FP(df_foldseek)
print('Foldseek 75% cov.:')
y_foldseek_cov, x_foldseek_cov = get_TP_FP(df_foldseek_cov)
print('Blastp all:')
y_blastp, x_blastp = get_TP_FP(df_blastp)


ax = sns.lineplot(x=x_blastp, y=y_blastp)
ax = sns.lineplot(x=x_foldseek, y=y_foldseek)
ax = sns.lineplot(x=x_foldseek_cov, y=y_foldseek_cov)

ax.legend(['RBH (75% cov.)', 'RBSH all', 'RBSH (75% cov.)'])
ax.set(xlabel = 'False positives', ylabel='True positives')
plt.savefig(args.png_out, dpi=600, bbox_inches='tight')
