import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--dom_results', required=True,
                    help='domain results files for RBSH hits')
parser.add_argument('--panther_p1', required=True,
                    help='panther classifications for organism of prot1 '
                    '(downloaded for protein ids of organism 1 from '
                    'http://pantherdb.org/geneListAnalysis.do)')
parser.add_argument('--panther_p2', required=True,
                    help='panther classifications for organism of prot2 '
                    '(downloaded for protein ids of organism 2 from '
                    'http://pantherdb.org/geneListAnalysis.do)')
parser.add_argument('--out_fh_name', required=True,
                    help='output file name for .csv and .xlsx files')
args = parser.parse_args()


def read_in_panther(fh, protein_id):
    cols = ['prot{}'.format(protein_id), 'Subfamily_p{}'.format(protein_id),
            'description_p{}'.format(protein_id)]
    df = pd.read_csv(fh, sep='\t', names=cols)
    df['Family_p{}'.format(protein_id)] = df['Subfamily_p{}'.format(protein_id)].apply(
        lambda x: x.split(':')[0])
    return df


df_prot1 = read_in_panther(args.panther_p1, 1)
df_prot2 = read_in_panther(args.panther_p2, 2)

df_dom = pd.read_csv(args.dom_results)
df_merged = df_dom.merge(
    df_prot1, on='prot1', how='left').merge(
        df_prot2, on='prot2', how='left')

df_merged['Panther'] = np.where(
    df_merged['Family_p1'] == df_merged['Family_p2'], 'Same', 'Diff')
df_merged.to_csv('{}.csv'.format(args.out_fh_name), index=False)
df_merged.to_excel('{}.xlsx'.format(args.out_fh_name), index=False)
