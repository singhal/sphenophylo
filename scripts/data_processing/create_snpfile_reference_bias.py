import re
import gzip
import argparse

parser = argparse.ArgumentParser(description = "make compact SNP output", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--lineage", type = str, default = None, help = "lineage for which to run script")
parser.add_argument("--dp", type = int, default = 5, help = "keep sites with >= this depth")
args = parser.parse_args()

lfile = '/media/babs/brains3/spheno/ref_genome/skink_SOAPdenovo_k107.scafSeq_sortedbylength.lengths.csv'
DEPTH = int(args.dp)
sp = args.lineage

f = open(lfile, 'r')
head = f.readline()
seq = {}
for l in f:
    d = re.split(',', l.rstrip())
    if int(d[1]) >= 5000:
        seq[d[0]] = int(d[1])
f.close()

def get_geno(geno):
    genonums = {'A': {'A': 1, 'T': 2, 'C': 3, 'G': 4},
                'T': {'A': 2, 'T': 5, 'C': 6, 'G': 7},
                'C': {'A': 3, 'T': 6, 'C': 8, 'G': 9},
                'G': {'A': 4, 'T': 7, 'C': 9, 'G': 0}}
    return genonums[geno[0]][geno[1]]
    
vcf = '/media/babs/brains3/spheno/ddRAD_genome/variants/%s.vcf.gz' % sp
f = gzip.open(vcf, 'r')
for l in f:
    l = l.decode("utf-8")
    if re.search('#CHROM', l):
        d = re.split('\t', l.rstrip())
        inds = d[9:len(d)]

        s = {}
        for id, slen in seq.items():
            s[id] = [[-9] * len(inds) for i in range(slen)]
            
    elif not re.search('#', l):
        d = re.split('\t', l.rstrip())
        c = d[0]
        if c in s:
            pos = int(d[1]) - 1
            a = [d[3]] + re.split(',', d[4])
            genos = d[9:len(d)]
            dpix = re.split(":", d[8]).index("DP")
            for ix, geno in enumerate(genos):
                geno = re.split(":", geno)
                dp = int(geno[dpix])
                if dp >= DEPTH:
                    geno = [a[int(x)] for x in re.split('/', geno[0])]
                    genonum = get_geno(geno)
                    s[c][pos][ix] = genonum
f.close()

out = re.sub('.vcf.gz', '.dp%s.snp' % DEPTH, vcf)
o = open(out, 'w')
o.write('# %s\n' % ' '.join(inds))
for c in s:
    for ix, pos in enumerate(s[c]):
        g = ['.' if x == -9 else str(x) for x in pos ]
        g = ''.join(g)
        if not re.search('^\.+$', g):
            o.write('%s,%s,%s\n' % (c, ix + 1, g))
o.close()
