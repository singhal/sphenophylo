import re
import subprocess
import os
import argparse
import random
import string
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Get homology file for groups of interest.')
parser.add_argument('--l', help="cluster column")
parser.add_argument('--n', help="cluster name")
parser.add_argument('--c', help="cluster percentage")
parser.add_argument('--d', help="drop num")
args = parser.parse_args()
name = args.n
col = args.l
dump = int(args.d)
c = float(args.c)
genus = name

def get_seq(seqfile):
	seq = {}
	f = open(seqfile, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()
	return seq

indfile = '/home/babs/spheno/sphenomorphine_ddRAD_v3.csv'
d = pd.read_csv(indfile)
d = d.loc[d.RECENT_RUN == False, ]
d = d.loc[ d[col] == name, ]

outdir = '/home/babs/spheno/phylogeny/ddRAD/homology/'
fadir = '/home/babs/spheno/phylogeny/ddRAD/ind_assemblies_derep/'
orgfiles = [(ind, '%s%s.fa' % (fadir, ind)) for ind in d.SAMPLE_ID.tolist()]
orgfiles = dict(orgfiles)

# cluster seqs at this value
WCLUST = c
# min number of seqs needed to include
minseq = 100

# dump percentage (locus needs to be in at least dumpper * dump individuals)
dumpper = 0.2

out = outdir + '%s_inds.txt' % name
o = open(out, 'w')
# filter files
files = {}
cts = {}
for ind, file in orgfiles.items():
	if os.path.isfile(file):
		seq = get_seq(file)
		numseq = len(seq)
		if numseq >= minseq:
			files[ind] = file
			cts[ind] = numseq
			o.write('%s\tKEEP\n' % ind)
		else:
			o.write('%s\tdropped_too_few_seq\n' % ind)	
	else:	
		o.write('%s\tdropped_no_FASTA\n' % ind)
o.close()

def create_starter(dir, file, genus, ix, ind):
	homhash = {}

	starting = '%s%s.tmp.fa' % (dir, genus)
	f = open(file, 'r')
	o = open(starting, 'w')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			newid = '%s_%s' % (genus, ix)

			homhash[newid] = {}
			homhash[newid][id] = '+'
			
			seq = f.readline().rstrip()
			o.write('>%s\n%s\n' % (newid, seq))
			ix += 1
	f.close()
	o.close()

	return starting, homhash, ix


def vsearch(dir, tmpfile, file, genus, num):
	out = '%s%s_%s_search' % (dir, genus, num)
	subprocess.call("vsearch --usearch_global %s --db %s --userout %s --id %s --userfields query+target+evalue+id+qstrand --strand both --threads 4" % (file, tmpfile, out, WCLUST), shell=True)

	return out


def create_new_tmp(dir, tmpfile, seqfile, results, homhash, genus, ix, ind, todump, dumpind):
	matches1 = {}
	matches2 = {}

	f = open(results, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		# is this step necessary?
		# makes sure it is 1 to 1
		match = d[0]
		if d[1] not in matches1 and match not in matches2:
			matches1[d[1]] = {'match': match, 'perc': float(d[3]), 'strand': d[4]}
			matches2[match] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
		elif match in matches2 and d[1] not in matches1:
			if float(d[3]) > matches2[match]['perc']:
				matches1[d[1]] = {'match': match, 'perc': float(d[3]), 'strand': d[4]}
				matches2[match] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
	f.close()
	os.remove(results)

	for c in matches2:
		homhash[matches2[c]['match']][c] = matches2[c]['strand']

	if todump:
		homhash2 = {}
		for c in homhash:
			nummatch = len(homhash[c])
			if nummatch >= dumpind:
				homhash2[c] = {}
				homhash2[c] = homhash[c]
	else:
		homhash2 = homhash
			

	# need to read in tmp file
	# only want to output into a new tmp file those things in homhash2
	# don't need to worry about the singletons in the ind file that did not get matched
	if todump:
		seq = get_seq(tmpfile)
		o = open(tmpfile, 'w')
		for c in homhash2:
			o.write('>%s\n%s\n' % (c, seq[c]))
		o.close()
	else:
		f = open(seqfile, 'r')
		o = open(tmpfile, 'a')
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l.rstrip()).group(1)
				seq = f.readline().rstrip()
				if id not in matches2:
					new_id = '%s_%s' % (genus, ix)
					ix += 1

					homhash2[new_id] = {}	
					homhash2[new_id][id] = '+'
			
					o.write('>%s\n%s\n' % (new_id, seq))
		f.close()
		o.close()

	return (tmpfile, homhash2, ix)

inds = sorted(cts, key=cts.get, reverse=True)
ix = 0
tmpfile, homhash, ix = create_starter(fadir, files[inds[0]], genus, ix, inds[0])
for num, ind in enumerate(inds[1:]):
	file = files[inds[num + 1]]
	results = vsearch(fadir, tmpfile, file, genus, num)
	print("we are on ind %s: %s (%s)" % (num, ind, cts[ind]))

	num2 = num + 1
	if num2 % dump == 0:
		todump = True
	else:
		todump = False
	dumpind = dumpper * num2
	(tmpfile, homhash, ix) = create_new_tmp(fadir, tmpfile, file, results, homhash, genus, ix, ind, todump, dumpind)
os.remove(tmpfile)

o = open('%s%s_homology_across_species.c%s.txt' % (outdir, genus, WCLUST), 'w')
o.write('contig\tmatches\tnumMatches\n')
for c, matches in homhash.items():
	matches = ['%s:%s' % (match, homhash[c][match]) for match in matches]
	o.write('%s\t%s\t%s\n' % (c, ','.join(matches), len(matches)))
o.close()
