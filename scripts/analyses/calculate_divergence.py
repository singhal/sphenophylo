import re
import argparse
import os
import gzip
import pandas as pd
import numpy as np

def fst_reich(counts, sample_sizes):
    counts1 = np.array(counts)[0]
    counts2 = np.array(counts)[1]

    sample_sizes1 = np.array(sample_sizes).astype('float')[0]
    sample_sizes2 = np.array(sample_sizes).astype('float')[1]

    h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
    h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))

    N = []
    D = []

    for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
        n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
        N.append(n)
        d = n + _h1 + _h2
        D.append(d)

    F = np.sum(N) / np.sum(D)

    return F

def calculate_fst(snps, inds, inds1, inds2):

    counts = {'sp1': [], 'sp2': []}
    sizes = {'sp1': [], 'sp2': []}

    for c in snps:
        for pos in snps[c]:
            geno = snps[c][pos]
            geno = dict(zip(inds, geno))

            a1 = [geno[ind] for ind in inds1]
            a1 = [bp for geno in a1 for bp in geno]
            a1 = [x for x in a1 if x != 'N']

            a2 = [geno[ind] for ind in inds2]
            a2 = [bp for geno in a2 for bp in geno]
            a2 = [x for x in a2 if x != 'N']

            # only want to work with positions where
            # there are two snps and two snps only
            if len(set(a1 + a2)) == 2:
                # this is the variable site
                # determined arbitarly
                # fst doesn't have minor / major (i think?)
                allele = list(set(a1 + a2))[0]

                if len(a1) == 0:
                    count1 = np.nan
                else:
                    count1 = a1.count(allele)
                counts['sp1'].append(count1)
                sizes['sp1'].append(len(a1))

                if len(a2) == 0:
                    count2 = np.nan
                else:
                    count2 = a2.count(allele)
                counts['sp2'].append(count2)
                sizes['sp2'].append(len(a2))

    # confirm sample sizes
    alleles = np.array([counts['sp1'], counts['sp2']])
    sizes = np.array([sizes['sp1'], sizes['sp2']])
    to_mask = np.any(np.isnan(alleles), axis=0)
    alleles = alleles[:, ~to_mask]
    sizes = sizes[:, ~to_mask]

    if len(alleles[0]) > 0:
        fst = fst_reich(alleles, sizes)
        counts = len(alleles[0])
    else:
        fst = np.nan
        counts = 0
            
    return fst, counts

def calc_div_prop(a1, a2):
    alleles = a1 + a2
    if len(set(alleles)) > 1:
        n_c = 0
        n_diff = 0
        for x in a1:
            for y in a2:
                n_c += 1
                if x != y:
                    n_diff += 1
        div = n_diff / float(n_c)
    else:
        div = 0

    return div

def calc_pi_prop(alleles):
    if len(set(alleles)) > 1:
        n_c = 0
        n_diff = 0
        for i in alleles:
            for j in alleles:
                n_c += 1
                if i != j:
                    n_diff += 1
        pi_prop = n_diff / float(n_c)
    else:
        pi_prop = 0

    return pi_prop

def calculate_dxy(snps, inds, inds1, inds2):
    dxy = {'pi1': 0, 'pi2': 0, 'diff': 0, 'denom': 0}

    for c in snps:
        for pos in snps[c]:
            geno = snps[c][pos]
            geno = dict(zip(inds, geno))

            a1 = [geno[ind] for ind in inds1]
            a1 = [bp for geno in a1 for bp in geno]
            a1 = [x for x in a1 if x != 'N']

            a2 = [geno[ind] for ind in inds2]
            a2 = [bp for geno in a2 for bp in geno]
            a2 = [x for x in a2 if x != 'N']

            # print(a1)
            # print(a2)
            # print('****')

	    # need at least two chromosomes
            if len(a1) > 1 and len(a2) > 1:
		# don't mess with multiallelics
                if len(set(a1 + a2)) < 3:
                    dxy['denom'] += 1
                    dxy['diff'] += calc_div_prop(a1, a2)
                    dxy['pi1'] += calc_pi_prop(a1)
                    dxy['pi2'] += calc_pi_prop(a2)
    return dxy

def get_seq(seqfile):
    f = open(seqfile, 'r')
    id = ''
    seq = {}
    for l in f:
        if re.search('>', l):
            id = re.search('>(\S+)', l.rstrip()).group(1)
            seq[id] = ''
        else:
            seq[id] += l.rstrip()
    f.close()
    return seq

def calc_seq_div(s1, s2):
    s1 = s1.upper()
    s2 = s2.upper()

    bases = ['A', 'C', 'G', 'T']

    diff = 0
    length = 0
    
    for a, b in zip(s1, s2):
        if a in bases and b in bases:
            length += 1
            if a != b:
                diff += 1

    if length > 0:
        div = diff / float(length)
    else:
        div = -1

    return div
    
def seq_divergence(seq1, seq2):
    n_c = 0
    diff = 0

    for x, s1 in seq1.items():
        for y, s2 in seq2.items():
            div = calc_seq_div(s1, s2)
            if div > 0:
                diff += div
                n_c += 1

    if n_c > 0:
        div = diff / float( n_c )
    else:
        div = np.nan

    return div

def calc_mt_div(seq, inds1, inds2):
    seq1 = {}
    for ind in inds1:
        if ind in seq:
            seq1[ind] = seq[ind]
    seq2 = {}
    for ind in inds2:
        if ind in seq:
            seq2[ind] = seq[ind]

    res = {}
    res['mt_within1'] = seq_divergence(seq1, seq1)
    res['mt_within2'] = seq_divergence(seq2, seq2)
    res['mt_btn'] = seq_divergence(seq1, seq2)

    return res 

    
parser = argparse.ArgumentParser(description = "calculate divergence between groups",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# then get group 1 & 2
parser.add_argument("--g1", type = str, default = None, help = "group 1 to compare")
parser.add_argument("--g2", type = str, default = None, help = "group 2 to compare")
parser.add_argument("--sp", type = str, default = None, help = "species of groups")
parser.add_argument("--file", type = str, default = None, help = "file with ind to group assignemnts")
args = parser.parse_args()

# all the metadata
d = pd.read_csv(args.file)

m = '../Sphenomorphine_speciation/phylogeny/mtDNA/sphenomorphine_cytb.aligned.edited.fasta'

# get inds in group 1 & 2
g1 = args.g1
g2 = args.g2
inds1 = d.loc[d['pop'] == g1, "ind"].tolist()
inds2 = d.loc[d['pop'] == g2, "ind"].tolist()

inds = inds1 + inds2

s = {}
vcff = 'data/variants/%s.dp5.miss0.5.vcf.recode.vcf' % args.sp
f = open(vcff, 'r')
for ix, l in enumerate(f):
    # l = l.decode("utf-8")
    if re.search('^#CHROM', l):
        d = re.split('\t', l.rstrip())
        indix = d[9:]
    elif not re.search('^#', l):
        d = re.split('\t', l.rstrip())
        c = d[0]
        pos = int(d[1]) - 1

        if c not in s:
            s[c] = {}
        if pos not in s[c]:
            s[c][pos] = [['N', 'N']] * len(inds)

        # alleles
        snps = [d[3]] + re.split(',', d[4])
        a = {}
        for ix, snp in enumerate(snps):
            a[str(ix)] = snp.upper()
        a['.'] = 'N'
    
        gens = [re.search('^(\S\S\S)', x).group(1) for x in d[9:]]
        gens = [re.split('/', x) for x in gens]

        for ind in inds:
            if ind in indix:
                gen = gens[ indix.index(ind) ]
                gen = [a[x] for x in gen]
                s[c][pos][ inds.index(ind) ] = gen
f.close()


outf = 'data/divergence/%s.%s.%s.csv' % (args.sp, g1, g2)
o = open(outf, 'w')
o.write('metric,species1,species2,submetric,value\n')
# now calculate statistics
print("calculating dxy")
dxy = calculate_dxy(s, inds, inds1, inds2)
for key in dxy:
    o.write('dxy,%s,%s,%s,%s\n' % (g1, g2, key, dxy[key]))
print("calculating fst")
fst, fstc = calculate_fst(s, inds, inds1, inds2)
o.write('fst,%s,%s,%s,%s\n' % (g1, g2, 'fst', fst))
o.write('fst,%s,%s,%s,%s\n' % (g1, g2, 'fst_counts', fstc))


print("calculating mtdna")
seq = get_seq(m)
res = calc_mt_div(seq, inds1, inds2)
for key, value in res.items():
     o.write('mtdxy,%s,%s,%s,%s\n' % (g1, g2, key, value))
o.close()
