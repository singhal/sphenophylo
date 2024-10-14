import re
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='simplify individual assembly file')
parser.add_argument('--i', help="individual to sample")
args = parser.parse_args()
ind = args.i

WCLUST = 0.98

seqgz = '/home/babs/spheno/phylogeny/ddRAD/ind_assemblies/%s.fa.gz' % ind
seqfile = re.sub('.gz', '', seqgz)
subprocess.call("gunzip %s" % seqgz, shell = True)

indtmp1 = seqfile + '_1'
indtmp2 = seqfile + '_2'
indtmp3 = seqfile + '_3'

subprocess.call("vsearch --derep_fulllength %s --output %s --fasta_width 0 --strand both" %  (seqfile, indtmp1), shell=True)
subprocess.call("vsearch --sortbylength %s --output %s" % (indtmp1, indtmp2), shell=True)
subprocess.call("vsearch --cluster_smallmem %s --centroids %s --id %s --usersort --fasta_width 0 --strand both --minsl 0.5 --query_cov 0.7" % (indtmp2, indtmp3, WCLUST), shell=True)

os.remove(indtmp1)
os.remove(indtmp2)

def get_seq(seqfile):
    f = open(seqfile, 'r')
    seq = {}
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l).group(1)
            seq[id] = ''
        else:
            seq[id] += l.strip()
    f.close()
    return seq

seq = get_seq(indtmp3)
new = '/home/babs/spheno/phylogeny/ddRAD/ind_assemblies_derep/%s.fa' % ind
f = open(new, 'w')
ix = 1
for id, s in seq.items():
    f.write('>%s_%s\n%s\n' % (ind, ix, s))
    ix = ix + 1
f.close()

os.remove(indtmp3)
