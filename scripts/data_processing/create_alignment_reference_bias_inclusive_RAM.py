import re
import gzip
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = "make concat SNP output", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--level", type = str, default = None, help = "cluster level")
parser.add_argument("--name", type = str, default = None, help = "name of cluster")
parser.add_argument("--miss", type = float, default = None, help = "missingness level")
parser.add_argument("--outsp", type = str, default = None, help = "outgroup")
parser.add_argument("--mito", action='store_true', default = None, help = "add flag if you want to include mito data")
args = parser.parse_args()

cl = args.level
clname = args.name
MISS = float(args.miss)
out = re.split(',', args.outsp)

print("**************")
print("%s: %s" % (cl, clname))

cfile = '/media/babs/brains3/spheno/metadata/sphenomorphine_all_v6.csv'

def read_seq(seqfile):

    id = ''
    seq = {}
    o = open(seqfile, 'r')
    for l in o:
        l = l.rstrip()
        if re.search('>', l):
            id = re.search('>(\S+)', l).group(1)
            seq[id] = ''
        else:
            seq[id] += l
    o.close()

    return seq

########
# start with getting the ddRAD data
########
dd = pd.read_csv(cfile)

d1 = dd.loc[dd['RECENT_RUN'] == False, ]
d1 = d1.loc[d1['ddRAD'] == True, ]

sps = d1.loc[d1[cl] == clname, 'CURRENT_TAXON'].unique().tolist()
inds = d1.loc[d1.CURRENT_TAXON.isin(sps), 'SAMPLE_ID'].tolist() + out

# add in outgroup
sps = sps + d1.loc[d1.SAMPLE_ID.isin(out), 'CURRENT_TAXON'].tolist()

snp = {'1': 'A', '2': 'W', '3': 'M', '4': 'R', '5': 'T',
       '6': 'Y', '7': 'K', '8': 'C', '9': 'S', '0': 'G',
       '.': '-'}

inds2 = inds
s = {}
for sp in sps:
    print("doing %s ddrad" % sp)
    sfile = '/media/babs/brains3/spheno/ddRAD_genome/variants/%s.dp5.snp' % sp
    f = open(sfile, 'r')
    i = f.readline().rstrip()
    i = re.split(' ', re.sub('^# ', '', i))

    indixs = [999] * len(i)
    for indix, x in enumerate(i):
        if x in inds2:
            indixs[indix] = inds2.index(x)
    
    for l in f:
        d = re.split(',', l.rstrip())
        c = d[0]
        pos = int(d[1]) - 1

        if c not in s:
            s[c] = {}
        if pos not in s[c]:
            s[c][pos] = ['-'] * len(inds2)
        
        gens = [snp[x] for x in list(d[2])]
        for indix, gen in zip(indixs, gens):
            if indix < 999:
                s[c][pos][indix] = gen
    f.close()

keep = {}
for c in s:
    for pos in s[c]:
        bp = s[c][pos]
        if bp.count('-') / float(len(bp)) < MISS:
            if c not in keep:
                keep[c] = []
            keep[c].append(pos)


contigs = list(keep.keys())
inds3 = {}
for ind in inds2:
    seq = ''
    indix = inds2.index(ind)
    for c in contigs:
        bps = sorted(keep[c])
        for pos in keep[c]:
            seq += s[c][pos][indix]
    inds3[ind] = seq


################
# now sqcl
################

print("doing sqcl!")
sqclfile = '/media/babs/brains3/spheno/SqCL/concatenated/concat_ind0.05_loci0.7_all_n25.fasta'
sqcl = read_seq(sqclfile)

# get rid of the two bad individuals
del sqcl['CUMV_14452_Le_bipe']
del sqcl['UMMZ_244315_ct_quat']

# add in names of two samples that have diff names across datasets
sqcl['NA_CCM6337_Ct_deca'] = sqcl['CCM6337']
sqcl['NA_CCM5506_Ct_deca'] = sqcl['CCM5506']

d2 = dd.loc[dd['SQCL'] == True, ]
s_inds = d2.loc[d2[cl] == clname, 'SAMPLE_ID'].unique().tolist()

# add in outgroup
s_inds2 = s_inds + out

sqclinds = {}
for s_ind in s_inds2:
    seq = sqcl[s_ind]
    sqclinds[s_ind] = seq


#######################
# now mtDNA
#######################

def run_mito():
    print("doing mito!")
    mtfile = '/media/babs/brains3/spheno/mtDNA/sphenomorphine_cytb.aligned.edited.fasta'
    mtdna = read_seq(mtfile)

    d3 = dd.loc[dd['CYTB'] == True, ]
    m_inds = d3.loc[d3[cl] == clname, 'SAMPLE_ID'].unique().tolist()

    # add in outgroup
    m_inds2 = m_inds + out
    mtinds = {}

    for m_ind in m_inds2:
        mtinds[m_ind] = mtdna[m_ind]

    return mtinds

if args.mito:
    mtinds = run_mito()
else:
    mtinds = {}
    
## now print it out!
if args.mito:
    outf = '/media/babs/brains3/spheno/ddRAD_genome/alignments/%s.%s.miss%s.all.fasta' % (clname, cl, MISS)
    partf = '/media/babs/brains3/spheno/ddRAD_genome/alignments/%s.%s.miss%s.all.partitions' % (clname, cl, MISS)
else:
    outf = '/media/babs/brains3/spheno/ddRAD_genome/alignments/%s.%s.miss%s.nuc.fasta' % (clname, cl, MISS)
    partf = '/media/babs/brains3/spheno/ddRAD_genome/alignments/%s.%s.miss%s.nuc.partitions' % (clname, cl, MISS)
    
newseq = {}
allinds = list(set(list(mtinds.keys()) + list( sqclinds.keys() ) + list(inds3.keys()))) 
seqlen = len(list(inds3.values())[0])
if len(mtinds) > 0:
    mtlen = len(list(mtinds.values())[0])
else:
    mtlen = 0
sqlen = len(list(sqclinds.values())[0])

print(seqlen, mtlen, sqlen)


p = open(partf, 'w')
p.write("DNA, part_ddrad = 1-%s\n" % (seqlen))
second = seqlen + sqlen
p.write("DNA, part_sqcl = %s-%s\n" % (seqlen + 1, second))
if len(mtinds) > 0:
    third = second + mtlen
    p.write("DNA, part_mtdna = %s-%s\n" % (second + 1, third))
p.close()

o = open(outf, 'w')
for ind in allinds:
    if ind in inds3:
        seq1 = inds3[ind]
    else:
        seq1 = '-' * seqlen

    if ind in sqclinds:
        seq2 = sqclinds[ind]
    else:
        seq2 = '-' * sqlen
        
    if ind in mtinds:
        seq3 = mtinds[ind]
    else:
        seq3 = '-' * mtlen

    seq = seq1 + seq2 + seq3

    newseq[ind] = seq
        
    if ind in out:
        ind = ind + '.outgroup'
    keep = False

    miss1 = seq1.count('-') / float(len(seq1))
    miss2 = len(seq2) - seq2.count('-')
    miss3 = len(seq3) - seq3.count('-')
                
    if miss3 > 500 or miss2 > 20000 or miss1 < 0.95:
        o.write('>%s\n%s\n' % (ind, seq))
    else:
        missper = seq.count('-') / float(len(seq))
        print("%s dropped due to high missing data (%.2f)" % (ind, missper))
o.close()

seqlen = len(list(newseq.values())[0])

print("%s: %s" % (len(inds3), seqlen))
print("**************")
