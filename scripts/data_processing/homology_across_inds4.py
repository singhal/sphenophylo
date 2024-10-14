import argparse
import glob
import os
import re

dir = '/home/babs/spheno/phylogeny/ddRAD/alignments/'
outdir = re.sub('alignments', 'concatenated', dir)

parser = argparse.ArgumentParser(description='make a concatenated file.')
parser.add_argument('--n', help="cluster name")
parser.add_argument('--m', help="missingness for a locus")
args = parser.parse_args()

subdir = os.path.join(dir, args.n)
alns = glob.glob(subdir + '/*aln')

def get_seq(seqfile):
	f = open(seqfile, 'r')
	id = ''
	indseq = {}
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			id = re.sub('_R_', '', id)
			indseq[id] = ''
		else:
			indseq[id] += l.rstrip()
	f.close()
	return indseq


def trim_seq(seqaln):
	seq = [True] * len(list(seqaln.values())[0])
	for ix, pos in enumerate(seq):
		sites = [s[ix] for c, s in seqaln.items()]
		permiss = sites.count('-') + sites.count('N') + sites.count('?') + sites.count('n')
		permiss = permiss / float(len(sites))
		if permiss > 0.7:
			seq[ix] = False
	seqaln2 = {}
	for id, s in seqaln.items():
		newseq = ''.join([bp if keep else '' for bp, keep in zip(s, seq)])
		seqaln2[id] = newseq.upper()
	return seqaln2

inds = {}
allseq = {}
for aln in alns:
	seq = get_seq(aln)
	seq2 = trim_seq(seq)
	for ind in seq2:
		if ind not in inds:
			inds[ind] = 1
	allseq[aln] = seq2

keepaln = []
for aln in allseq:
        numind = len(allseq[aln])
        perind = float(numind) / len(inds)
        if perind >= float(args.m):
                keepaln.append(aln)

mfile = os.path.join(outdir, '%s.concatenated%s.csv' % (args.n, args.m))
m = open(mfile, 'w')
m.write('group,individual,num_loci,percent_missing,drop\n')
outfile = os.path.join(outdir, '%s.concatenated%s.fa' % (args.n, args.m))
o = open(outfile, 'w')
for ind in inds:
	indseq = ''
	nloci = 0
	for aln in keepaln:
		if ind in allseq[aln]:
			indseq += allseq[aln][ind]
			nloci += 1
		else:
			seqlen = len(list(allseq[aln].values())[0])
			indseq += '-' * seqlen
	miss = indseq.count('-') / float(len(indseq))
	if miss > 0.95:
		print("dropping %s because of excess missing" % ind)
		drop = True
	else:
		o.write('>%s\n%s\n' % (ind, indseq))
		drop = False
	m.write('%s,%s,%s,%s,%s\n' % (args.n, ind, nloci, miss, drop))
o.close()
m.close()
