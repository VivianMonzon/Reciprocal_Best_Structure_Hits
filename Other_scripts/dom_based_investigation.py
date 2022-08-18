import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dom_tbl_prot1', required=True,
                    help='hmmsearch Pfam tbl out for prot1')
parser.add_argument('--dom_tbl_prot2', required=True,
                    help='hmmsearch Pfam tbl out for prot2')
parser.add_argument('--recis_foldseek', required=True,
                    help='recis ids only found with RBSH approach')
parser.add_argument('--fh_out_name', required=True,
                    help='output file name for csv and xlsx')
args = parser.parse_args()

cols = ['UniProtID', 'acc', 'len', 'query', 'acc_q', 'qlen', 'E-value',
        'score', 'bias', '#', 'of', 'c-Evalue', 'i-Evalue', 'score_dom',
        'bias_dom', 'from_a', 'to_a', 'from_b', 'to_b', 'from', 'to',
        'x', 'descript']


def dom_per_proteome(hmm_tbl):
    df = pd.read_csv(hmm_tbl, comment='#', sep=r"\s+", names=cols)
    df = df[['UniProtID', 'len', 'query', 'acc_q', 'from', 'to']]
    df['acc_q'] = df.acc_q.apply(lambda x: x.split('.')[0])
    doms_per_protein = df[['UniProtID', 'acc_q']].drop_duplicates()
    doms_per_proteome = doms_per_protein['acc_q'].value_counts().rename_axis(
        'domain').reset_index(name='counts')
    doms_per_proteome = dict(zip(
        doms_per_proteome.domain, doms_per_proteome.counts))
    return df, doms_per_proteome


def get_doms_per_protein(df):
    ids = list(set(df.UniProtID.tolist()))
    doms_per_protein = {}
    for	p in ids:
        df_one = df[df.UniProtID == p]
        doms_per_protein[p] = list(set(df_one['acc_q'].tolist()))
    return doms_per_protein, ids


def combine_doms_per_protein(dom_df, protein_ids, dom_per_proteome):
    df_out = pd.DataFrame(columns={'UniProtID', 'Domains', 'Doms_number'})
    for p in protein_ids:
        df_one = dom_df[dom_df.UniProtID == p]
        df_one['acc_q'] = df_one['acc_q'].apply(lambda x: x + ' (' + str(
            dom_per_proteome[x]) + ')')
        doms_per_protein = list(set(df_one.acc_q.tolist()))
        doms_per_protein_str = ','.join(doms_per_protein)
        results_one = {'UniProtID': p, 'Domains': doms_per_protein_str,
                       'Doms_number': len(doms_per_protein)}
        df_out = df_out.append(results_one, ignore_index=True)
    return df_out
        
        
dom_df_prot2, doms_per_proteome_prot2 = dom_per_proteome(args.dom_tbl_prot2)
dom_df_prot1, doms_per_proteome_prot1 = dom_per_proteome(args.dom_tbl_prot1)


doms_per_protein_prot2, prot2_p_w_doms = get_doms_per_protein(dom_df_prot2)
doms_per_protein_prot1, prot1_p_w_doms = get_doms_per_protein(dom_df_prot1)

doms_adapt_prot2 = combine_doms_per_protein(dom_df_prot2, prot2_p_w_doms,
                                            doms_per_proteome_prot2)
doms_adapt_prot1 = combine_doms_per_protein(dom_df_prot1, prot1_p_w_doms,
                                            doms_per_proteome_prot1)

df_recis = pd.read_csv(args.recis_foldseek, names=['Recis'])
df_recis['prot1'] = df_recis['Recis'].apply(lambda x: x.split('_')[0])
df_recis['prot2'] = df_recis['Recis'].apply(lambda x: x.split('_')[1])
recis_dict = dict(zip(df_recis.prot1, df_recis.prot2))
all_prot1 = df_recis['prot1'].tolist()
all_prot2 = df_recis['prot2'].tolist()
prot1_wo_doms = [x for x in all_prot1 if x not in prot1_p_w_doms]
prot2_wo_doms = [x for x in all_prot2 if x not in prot2_p_w_doms]
prot1_wo_doms_dict = {}
for p1 in prot1_wo_doms:
    prot1_wo_doms_dict[p1] = ['P1']

prot2_wo_doms_dict = {}
for p2 in prot2_wo_doms:
    prot2_wo_doms_dict[p2] = ['P2']

doms_per_protein_prot2_all = {**doms_per_protein_prot2, **prot2_wo_doms_dict}
doms_per_protein_prot1_all = {**doms_per_protein_prot1, **prot1_wo_doms_dict}

recis_common_doms = {}
for k, v in recis_dict.items():
    common_doms = len([x for x in doms_per_protein_prot1_all[k] if x in
                       doms_per_protein_prot2_all[v]])
    recis = k + '_' + v
    recis_common_doms[recis] = common_doms
df_recis_common_doms = pd.DataFrame(recis_common_doms.items(),
                                    columns=['Recis', 'Same_dom'])

df_recis = df_recis.merge(df_recis_common_doms, on='Recis', how='left')

df_recis_dom_info = df_recis.merge(
    doms_adapt_prot1, left_on='prot1', right_on='UniProtID', how='left').merge(
        doms_adapt_prot2, left_on='prot2', right_on='UniProtID', how='left')

df_recis_dom_info = df_recis_dom_info.rename(
    {'UniProtID_x': 'UniProtID_p1', 'Doms_number_x': 'Doms_number_p1',
     'Domains_x': 'Domains_p1', 'UniProtID_y': 'UniProtID_p2',
     'Doms_number_y': 'Doms_number_p2', 'Domains_y': 'Domains_p2'}, axis=1)
df_recis_dom_info.to_csv('{}.csv'.format(args.fh_out_name), index=False)
df_recis_dom_info.to_excel('{}.xlsx'.format(args.fh_out_name), index=False)
