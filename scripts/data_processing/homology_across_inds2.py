import re
import subprocess
import os
import argparse
import random
import string
import sys
import pandas as pd
import glob

parser = argparse.ArgumentParser(description='get outgroup for homology file')
parser.add_argument('--l', help="cluster column")
parser.add_argument('--i', help="cluster ingroup name")
parser.add_argument('--o', help="cluster outgroup name")
args = parser.parse_args()

ingen = args.i
outgen = args.o
dir = '/home/babs/spheno/phylogeny/ddRAD/'
dfile = '/home/babs/spheno/sphenomorphine_ddRAD_v3.csv'
d = pd.read_csv(dfile)
d = d.loc[d.RECENT_RUN == False, ]

outcol = {'CURRENT_TAXON': 'CLUSTER_SPECIESCOMPLEX', 'CLUSTER_SPECIESCOMPLEX': 'CLUSTER_SUBGENUS', 'CLUSTER_SUBGENUS': 'CLUSTER_GENUS', 'CLUSTER_GENUS': 'CLUSTER_GENUS'}

def get_inds(genus):
    indfile = os.path.join(dir, 'homology', '%s_inds.txt' % genus)
    inds = {}
    o = open(indfile, 'r')
    for l in o:
        d = re.split('\t', l.rstrip())
        if d[1] == "KEEP":
            inds[d[0]] = 1
    o.close()
    return inds

def get_seq(indseq, seqfile):
    f = open(seqfile, 'r')
    id = ''
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l).group(1)
            indseq[id] = ''
        else:
            indseq[id] += l.rstrip()
    f.close()
    return indseq

def ind_seq(inds):
    seq = {}
    for ind in inds:
        seqfile =  os.path.join(dir, 'ind_assemblies_derep', '%s.fa' % ind)
        seq = get_seq(seq, seqfile)
    return seq


def num_match(inds):
    if len(inds) == 2:
        num_match = 2
    elif len(inds) == 3:
        num_match = 2
    else:
        num_match = len(inds) *	0.5
    return num_match

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement
    
def get_match(genus, num, seq):
    files = glob.glob(dir + 'homology/' + '*homology_across_species*txt')
    genusmatch = '%s_homology' % genus
    indfile = [file for file in files if re.search(genusmatch, file)][0]

    all = {}
    o = open(indfile, 'r')
    head = o.readline()
    for l in o:
        d = re.split('\t', l.rstrip())
        if int(d[2]) >= num:
            m = re.split(',', d[1])
            orrs = [re.search(':(\S)', x).group(1) for x in m]
            inds = [re.search('(\S+):', x).group(1) for x in m]

            seqlen = dict([(ind, len(seq[ind])) for ind in inds])
             
            ind = max(seqlen, key=seqlen.get)
            ix = inds.index(ind)
            orr = orrs[ix]

            indseq = seq[ind].upper()
            if orr == '-':
                # reverse complement
                indseq = rev_comp(indseq)
            all[d[0]] = indseq
    o.close()
    return all, indfile
            

def print_seq(gen, seq):
    out = os.path.join(dir, 'homology', '%s_cluster.fa' % gen)
    o = open(out, 'w')
    for id, s in seq.items():
        o.write('>%s\n%s\n' % (id, s))
    o.close()
    return out

def run_blat(inf, out):
    blatout = '%s_%s.vsearch.out' % (ingen, outgen)
    subprocess.call("vsearch --usearch_global %s --db %s --userout %s --id 0.8 --userfields query+target+evalue+id+qstrand --strand both --threads 4" % (inf, out, blatout), shell = True)

    m = {}
    o = open(blatout, 'r')
    for l in o:
        d = re.split('\t', l.rstrip())
        if d[0] not in m:
            m[d[0]] = {'m': d[1], 'orr': d[4]}

    o.close()
    os.remove(blatout)
    return m

def print_aln(genus, num, seq, m, outseq, indfile):
    outdir = os.path.join(dir, 'alignments', '%s' % genus)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    o = open(indfile, 'r')
    head = o.readline()
    for l in o:
        d = re.split('\t', l.rstrip())
        if int(d[2]) >= num:
            mm = re.split(',', d[1])
            orrs = [re.search(':(\S)', x).group(1) for x in mm]
            inds = [re.search('(\S+):', x).group(1) for x in mm]

            sf = os.path.join(outdir, '%s.fa' % d[0])
            out = open(sf, 'w')
            for ix, ind in enumerate(inds):
                indseq = seq[ind].upper()
                if orrs[ix] == '-':
                    indseq = rev_comp(indseq)
                ind2 = re.sub('_\d+$', '', ind) 
                out.write('>%s\n%s\n' % (ind2, indseq))
            # print outgroup
            for outind in m:
                if d[0] in m[outind]:
                    # get seq
                    outid = m[outind][d[0]]['m']
                    seqout = outseq[outind][outid]
                    if m[outind][d[0]]['orr'] == '-':
                        seqout = rev_comp(seqout)
                    out.write('>%s\n%s\n' % (outind, seqout))
            out.close()
            
    
ininds = get_inds(ingen)
outinds = d.loc[ d[outcol[args.l]] == outgen, ]
outinds = outinds.loc[ outinds[args.l] != ingen, 'SAMPLE_ID'].tolist()

# do 5 outgroup for sequence
cts = {}
for outind in outinds:
    seq = {}
    outseq = os.path.join(dir, 'ind_assemblies_derep', '%s.fa' % outind)
    seq = get_seq(seq, outseq)
    if len(seq) > 20000:
        cts[outind] = 1
outinds = random.sample(sorted(cts), 5)

# get a dict of all the in-group inds with their seqs
inseq = ind_seq(ininds)

# number contigs are we retaining (50% or more)
innum = num_match(ininds)

# which contigs are these?
inall, indfile = get_match(ingen, innum, inseq)

# print out seqfiles for in
inout = print_seq(ingen, inall)

# run blat (hahah vsearch)
m = {}
outseq = {}
for outind in outinds:
    outseqfile = os.path.join(dir, 'ind_assemblies_derep', '%s.fa' % outind)
    m[outind] = run_blat(inout, outseqfile)
    tmpseq = {}
    tmpseq = get_seq(tmpseq, outseqfile)
    outseq[outind] = tmpseq

# add outgroups
# ideal output is one file per locus incl outgroup seq
print_aln(ingen, innum, inseq, m, outseq, indfile)
