import pandas as pd
from urllib import request
import gzip
from Bio.SeqUtils import IUPACData
from Bio import SeqIO
from io import StringIO

hist_file = str(snakemake.input['history'])
pdb_files = snakemake.params['pdb_files']
output_file = str(snakemake.output)

accession = str(snakemake.wildcards['acc'])
date_max = pd.to_datetime(snakemake.config['uniprot']['date_max'])
host = str(snakemake.config['unisave']['host'])

residues = { k.upper(): v for k, v in IUPACData.protein_letters_3to1.items() }
residues['XAA'] = 'X'

def parse_pdb(pdb_file):
    start = stop = -1
    chain_seq = ''
    confs = []
    with gzip.open(pdb_file, 'rt') as pdb_fh:
        for line in pdb_fh:
            if line.startswith('ATOM'):
                atom = line[13:17].strip()
                chain = line[20:22].strip()
                if atom == 'N' and chain == 'A':
                    res = line[17:20]
                    chain_seq += residues[res]
            elif start < 0 and line.startswith('DBREF'):
                values = line.split()
                chain = values[2]
                if chain == 'A':
                    start = int(values[8]) - 1
                    stop  = int(values[9])
    return chain_seq, start, stop

def get_ver(version):
    url = "{host}/{acc}?format={format}&versions={ver}".format(host = host, acc = accession, format = "txt", ver = version)
    with request.urlopen(url) as fp:
        txt = fp.read().decode('utf-8')
    return txt

seq_dict = {}
for pdb_file in pdb_files:
    pdb_seq, start, stop = parse_pdb(pdb_file)
    for i in range(stop - start):
        seq_dict[i + start] = pdb_seq[i]
seq = ''.join(v for k, v in sorted(seq_dict.items()))

data = pd.read_csv(hist_file, sep = "\t")
data['Date'] = pd.to_datetime(data['Date'])

versions = data[data['Date'] <= date_max].drop_duplicates('Sequence version')

found = False
for index, row in versions.iterrows():
    txt = get_ver(row['Entry version'])
    swiss = SeqIO.parse(StringIO(txt), 'swiss')
    record = next(swiss)
    if record.seq == seq:
        found = True
        break
assert found, "Matching sequence for pdb {acc} not found".format(acc = accession)
with open(output_file, 'w') as fp:
    fp.write(txt)
