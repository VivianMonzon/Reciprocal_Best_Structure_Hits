import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--cif_file', required=True,
                    help='cif files from downloaded with pdb files '
                    'from AlphaFold database')
parser.add_argument('--out_seq_fa', required=True,
                    help='Novel file containing the protein sequences')
args = parser.parse_args()

listings = []
lines = [line.rstrip('\n') for line in open(args.cif_file)]
ignore = False

for line in lines:
    if line.startswith('_entity_poly.pdbx_seq_one_letter_code'):
        ignore = False
    elif line.startswith('_entity_poly.pdbx_strand_id'):
        ignore = True
    elif line.startswith('_entry.id AF-'):
        ignore = True
    if not ignore:
        listings.append(line.strip())

protein_id = listings[0].split('-')[1].split('-')[0]

seq = ''
for li in listings[2:]:
    seq += li
    if li == ';':
        break

seq = seq.split('_entity_poly.pdbx_seq_one_letter_code', 1)[1]
if ';' in seq:
    seq = seq.split(';')[1].split(';')[0]
if '_' in seq:
    seq = seq.split('_')[0]
seq = seq.replace(' ', '')

fh_out = open(args.out_seq_fa, 'a')
fh_out.write('>{}\n{}\n'.format(protein_id, seq))
