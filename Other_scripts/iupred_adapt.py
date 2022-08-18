import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fh_in')
parser.add_argument('--fh_out')
args = parser.parse_args()


df_iupred = pd.read_csv(args.fh_in, sep='\t',
                        names=['SEQID', 'POS', 'RES', 'IUPRED2'], comment='#')
df_iupred['ID'] = df_iupred.SEQID.apply(
    lambda x: x.split('|')[1].split('|')[0] if '|' in x else x)
protein_ids = list(set(df_iupred.ID.tolist()))
max_dict = {}
mean_dict = {}
fraction_disorder = {}
for x in protein_ids:
    df_one = df_iupred[df_iupred.ID == x]
    df_one_len = df_one.shape[0]
    df_one_max = round(df_one.IUPRED2.max(), 2)
    df_one_mean = round(df_one.IUPRED2.mean(), 2)
    perc_dis = round(((df_one[df_one.IUPRED2 > 0.5].shape[0]) / (df_one_len)) * 100, 2)
    max_dict[x] = df_one_max
    mean_dict[x] = df_one_mean
    fraction_disorder[x] = perc_dis
iupred_max = pd.DataFrame(max_dict.items(),
                          columns=['ID', 'iupred_max'])
iupred_mean = pd.DataFrame(mean_dict.items(),
                           columns=['ID', 'iupred_mean'])
iupred_frac = pd.DataFrame(fraction_disorder.items(),
                           columns=['ID', 'frac_disordered'])
df_iupred = df_iupred.merge(
    iupred_max, how='outer').merge(
        iupred_mean, how='outer').merge(
            iupred_frac, how='outer')
df_iupred = df_iupred[['ID', 'iupred_max', 'iupred_mean', 'frac_disordered']]
df_iupred = df_iupred.drop_duplicates()
df_iupred.to_csv(args.fh_out, index=False)
