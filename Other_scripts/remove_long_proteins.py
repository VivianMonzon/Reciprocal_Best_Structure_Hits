import pandas as pd
import os

df_ids = pd.read_csv('structure_db/H_sapiens_AF_ids_w_long_proteins.csv',
                     names=['ID'])
ids = df_ids.ID.tolist()

long_proteins = list(set([x for x in ids if ids.count(x) > 1]))
fh_out = open('data/H_sapiens_long_proteins.csv', 'w')
for p in long_proteins:
    fh_out.write('{}\n'.format(p))
    
remain_proteins = [x for x in ids if x not in long_proteins]

for rp in remain_proteins:
    os.system(
        'cp structure_db/'
        'H_sapiens_UP000005640_9606_HUMAN_v2/AF-{}-F1-model_v2.*.gz '
        'structure_db/H_sapiens_wo_long_proteins/'.format(rp))

