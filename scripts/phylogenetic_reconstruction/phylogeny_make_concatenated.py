import re
import os
import glob
import argparse
import random
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indcutoff", required=True, help='percent of loci needed for ind to keep ind')
parser.add_argument("-l", "--loccutoff", required=True, help='percent of inds needed for loci to keep loci')
parser.add_argument("-a", "--alndir", required=True, help='directory with alignments')
parser.add_argument("-o", "--outdir", required=True, help="output directory for astral trees")
parser.add_argument("-n", "--number", required=False, help="number of loci desired; leave blank if all")
parser.add_argument("-d1", "--dropind", required=False, help="csv file with taxa to drop; first column has names")
parser.add_argument("-d2", "--droploci", required=False, help="csv file with loci to drop; first column has names")
parser.add_argument("-k", "--keeploci", required=False, help="csv file with loci to keep; first column has names")
parser.add_argument("--AHE", action="store_true", default=False, help='add flag if you only want AHE loci')
parser.add_argument("--uce", action="store_true", default=False, help='add flag if you only want uce loci')
args = parser.parse_args()

indcutoff = float(args.indcutoff)
loccutoff = float(args.loccutoff)

# filter inds with a lot of missing data
def inds_drop(inds, alns2, cutoff):
	to_drop = []
	for ind in inds:
		count = 0
		for locus in alns2:
			if ind in alns2[locus]:
				count += 1
		amt = count / float(len(alns2))
		if amt < cutoff:
			to_drop.append(ind)
	return to_drop

# filter loci that have a lot of missing data
def loci_drop(loci, inds, alns2, cutoff):
	to_drop = []
	for locus in loci:
		count = 0
		for ind in alns2[locus]:
			if ind in inds:
				count += 1
		amt = count / float(len(inds))
		if amt < cutoff:
			to_drop.append(locus)
	return to_drop

def get_inds(aln):
	inds = []
	f = open(aln, 'r')
	for l in f:
		if re.search('>', l):
			ind = re.search('>(\S+)', l).group(1)
			ind = re.sub('^_R_', '', ind)
			inds.append(ind)
	f.close()
	return inds

# first need to create a hash of individuals and loci
def get_data(alns):
	inds = {}
	loci = {}
	alns2 = {}
	
	for ix, aln in enumerate(alns):
		locus = re.search('([^/]+).fasta', aln).group(1)
		
		loci[locus] = 1
		tmpinds = get_inds(aln)
		for ind in tmpinds:
			if ind not in inds:
				inds[ind] = 1					
		alns2[locus] = tmpinds

	return alns2, inds, loci

def get_seq(aln):
	f = open(aln, 'r')
	tmpseq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			id = re.sub('^_R_', '', id)
			tmpseq[id] = ''
		else:
			tmpseq[id] += l.rstrip()
	f.close()

	loclen = len(list(tmpseq.values())[0])

	return tmpseq, loclen

seqtype = 'all'
alndir = args.alndir
if args.AHE:
	seqtype = 'AHE'
	alns = glob.glob(alndir + '/AHE*aln')
elif args.uce:
	seqtype = 'uce'
	alns = glob.glob(alndir + '/uce*aln')
else:
	alns = glob.glob(alndir + '/*aln')

# revise the inds and the loci with drop data
alns2, inds, loci = get_data(alns)
if args.dropind:
	xx = pd.read_csv(args.dropind)
	drop_inds = xx.loc[:, 0].tolist()
else:
	drop_inds = []
if args.droploci:
	xx = pd.read_csv(args.droploci)
	drop_loci = xx.loc[:, 0].tolist()
else:
	drop_loci = []

if args.keeploci:
        xx = pd.read_csv(args.keeploci)
        loci = xx.locus.tolist()
else:
        loci = [x for x in loci if x not in drop_loci]
inds = [x for x in inds if x not in drop_inds]

# identify the drop_inds
drop_inds = inds_drop(inds, alns2, indcutoff)
# revise the inds
inds = [x for x in inds if x not in drop_inds]
# revise the drop loci
drop_loci = loci_drop(loci, inds, alns2, loccutoff)
loci = [x for x in loci if x not in drop_loci]

allseq = {}
for x in inds:
	allseq[x] = ''
if args.number:
	loci = random.sample(loci, int(args.number))
alns2 = []
	
out = os.path.join(args.outdir, 'concat_ind%s_loci%s_%s_n%s.phy' % (indcutoff, loccutoff, seqtype, len(loci)))
partfile = os.path.join(args.outdir, 'concat_ind%s_loci%s_%s_n%s.partitions' % (indcutoff, loccutoff, seqtype, len(loci)))
o = open(out, 'w')

pout = open(partfile, 'w')
# partdir = re.sub('trim', 'partitions', alndir)
totlen = 0

def get_partition(partdir, locus, totlen, parts):
	pfile = os.path.join(partdir, '%s.partitions' % locus)
	f = open(pfile, 'r')
	for line in f:
		if re.search('uce', locus):
			if re.search('p1', line):
				key = 'nc'
			else:
				key = 'uce'
		else:
			key = re.search(', (\S+)=', line).group(1)
	
		p = re.search('=(\S+)', line).group(1)
		p = re.split(',', p)
		pnew = []
		for x in p:
			if re.search('\\\\3', x):
				x = re.sub('\\\\3', '', x)
				vals = [int(x) + totlen for x in re.split('-', x)]
				pnew.append('%s-%s\\3' % (vals[0], vals[1]))
			else:
				vals = [int(x) + totlen for x in re.split('-', x)]
				pnew.append('%s-%s' % (vals[0], vals[1]))
	
		if len(parts[key]) == 0:
			parts[key] = ','.join(pnew)
		else:
			parts[key] = parts[key] + ',' + ','.join(pnew)
		
	return parts

# go through each aln
for ix, locus in enumerate(loci):
	print(ix)
	alnfile = os.path.join(alndir, '%s.fasta.aln' % locus)
	seq, loclen = get_seq(alnfile)
	# parts = get_partition(partdir, locus, totlen, parts)
	pout.write('%s, %s=%s\n' % (locus, totlen, totlen + loclen - 1))
	totlen += loclen
	for ind in allseq:
		if ind in seq:
			allseq[ind] += seq[ind]
		else:
			allseq[ind] += '-' * loclen
pout.close()

print('number of loci: %s' % len(loci))
print('number of inds: %s' % len(allseq))
o.write('%s\t%s\n' % (len(allseq), totlen))
# make alignment
for ind in allseq:
	o.write('%s\t%s\n' % (ind, allseq[ind]))
o.close()
