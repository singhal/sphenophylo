import re
import os
import argparse
import random

parser = argparse.ArgumentParser(description="create pop gen files")
parser.add_argument('--structure',
			default=False,
			action='store_true',
			help="generate structure input")
parser.add_argument('--admixture',
			default=False,
			action='store_true',
			help="generate admixture input")
parser.add_argument('--adegenet',
			default=False,
			action='store_true',
			help="generate adegenet imnput")
parser.add_argument('--svd',
			default=False,
			action='store_true',
			help = "generate svd ind file")
parser.add_argument('--random',
			default=False,
			action='store_true',
			help="keep just one snp. default is to use all.")
parser.add_argument('--stem', 
			default=None,
			help="stem for all output files; can just default")
parser.add_argument('--MISS',
			default = None,
			help = "percent missing")
parser.add_argument('--MAC',
			default = None,
			help = "mininum allele count req'd; rec 2")
parser.add_argument('--sample',
			default = None,
			help = "# of snps to randomly sample")
parser.add_argument('--vcffile',
			default = None,
			help = "input vcf file")
parser.add_argument('--outdir',
			default = None,
			help = "output dir")
parser.add_argument('--drop',
			default = None,
			help = "individuals to drop")
args = parser.parse_args()

snps = args.vcffile
outdir = args.outdir
if not os.path.isdir(outdir):
	os.mkdir(outdir)

def make_code():
	ix = 0
	code = {}
	invcode = {}

	allbp = ['A', 'T', 'C', 'G']

	for i, bp1 in enumerate(allbp):
		if bp1 not in invcode:
			invcode[bp1] = {}
		for bp2 in allbp[i:]:
			if bp2 not in invcode:
				invcode[bp2] = {}
			code[str(ix)] = [bp1, bp2]
			invcode[bp1][bp2] = str(ix)
			invcode[bp2][bp1] = str(ix)
			ix += 1

	invcode['N'] = {}
	invcode['N']['N'] = '-'
	code['-'] = ['N', 'N']
	
	return code, invcode


def get_af(snp, code):
	a = [code[x] for x in snp]
	a = [x for geno in a for x in geno]
	a = [x for x in a if x != 'N']

	uniq = list(set(a))
	cts = [ a.count(bp) for bp in uniq]

	return len(uniq), min(cts)


def get_random(var):
	var2 = {}
	for c in var:
		keep = random.choice(list(var[c].keys()))
		var2[c] = {}
		var2[c][keep] = var[c][keep]

	return var2


def parse_snps(snps, invcode, MISS, MAC, drop):
	keep = {}

	f = open(snps, 'r')
	
	for l in f:
		if re.search('#CHROM', l):
			d = re.split('\t', l.rstrip())
			inds = d[9:]
		if not re.search('^#', l):
			d =	re.split('\t', l.rstrip())
			snps = d[9:]

			good = True
			# get rid of indels & non biallelics
			if len(d[3]) > 1 or len(d[4]) > 1:
				good = False

			snps2 = []
			for snp, ind in zip(snps, inds):
				if ind not in drop:
					snps2.append( re.search('(^\S\S\S)', snp).group(1) )
			complete = 1 - (snps2.count('./.') / float(len(snps2)))

			snps = [re.search('(^\S\S\S)', snp).group(1)  for snp in snps]

			if complete > MISS and good:
				a = [re.split('/', snp) for snp in snps]
				a = [x for b in a for x in b]
				af1 = a.count('0')
				af2 = a.count('1')

				if af1 >= MAC and af2 >= MAC:
					# keep[d[0]][d[1]] = snp
					if re.search('_(\d+)', d[0]):
						chr = re.search('_(\d+)', d[0]).group(1)
					else:
						chr = d[0]
					if chr not in keep:
						keep[chr] = {}

					snp = []
					alleles = {'0': d[3], '1': d[4], '.': 'N'}
					for geno, ind in zip(snps, inds):
						geno = [alleles[x] for x in re.split('/', geno)]
						snp.append( invcode[geno[0]][geno[1]] )
			
					keep[chr][int(d[1])] = ''.join(snp)
						
	f.close()

	return keep, inds

def make_svd(inds, var, stem, code, drop):
	out = os.path.join(outdir, stem + '.nex')
	o = open(out, 'w')

	snp = {	'A': {'A': 'A', 'T': 'W', 'C': 'M', 'G': 'R'},
			'T': {'A': 'W', 'T': 'T', 'C': 'Y', 'G': 'K'},
			'C': {'A': 'M', 'T': 'Y', 'C': 'C', 'G': 'S'},
			'G': {'A': 'R', 'T': 'K', 'C': 'S', 'G': 'G'},
			'N': {'N': '-'}}

	seq = {}
	for ind in inds:
		seq[ind] = ''
	
	for c in sorted(var.keys()):
		for bp in sorted(var[c].keys()):
			for ind, geno in zip(inds, list(var[c][bp])):
				seq[ind] += snp[code[geno][0]][code[geno][1]]

	for ind in drop:
		if ind in seq:
			del seq[ind]
	maxlen = max([len(ind) for ind in seq]) + 1
	
	o.write("#NEXUS\n")
	o.write("\n")
	o.write("BEGIN DATA;\n")
	o.write("\tDIMENSIONS NTAX=%s NCHAR=%s;\n" % (len(seq), len(list(seq.values())[0])))
	o.write("\tFORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
	o.write("\tMATRIX\n")

	for ind in seq:
		s = seq[ind]
		indout = re.sub('-', '_', ind) + ' ' * (maxlen - len(ind))
		o.write('%s %s\n' % (indout, s))
	o.write(";\nEND;\n")
	o.close()


def make_structure(inds, b, stem, code, drop):
	out = os.path.join(outdir, '%s.str' % stem)
	o = open(out, 'w')

	for ix, ind in enumerate(inds):
		gen1 = [ind, '1', '0', '1', '1']
		gen2 = [ind, '1', '0', '1', '1']

		for c in sorted(b.keys()):
			for pos in sorted(b[c].keys()):
				gen1.append(b[c][pos][ix][0])
				gen2.append(b[c][pos][ix][1])

		if ind not in drop:
			o.write(' '.join(gen1) + '\n')
			o.write(' '.join(gen2) + '\n')

	o.close()

def make_binary(var, code):
	binary = {}
	for c in var:
		binary[c] = {}
		for pos in var[c]:
			genos = list(var[c][pos])
			genos = [code[geno] for geno in genos]
		
			alleles = [bp for geno in genos for bp in geno]
			alleles = [bp for bp in alleles if bp != 'N']
			alleles = list(set(alleles))
		
			a = {}
			for ix, allele in enumerate(alleles):
				a[allele] = str(ix)
			a['N'] = '-9'

			new = []
			for geno in genos:
				new.append([a[geno[0]], a[geno[1]]])
			binary[c][pos] = new
	return binary

def remove_miss(var, inds):
	tot_len = 0
	indmiss = dict([(ind, 0) for ind in inds])

	for c in var:
		for pos in var[c]:
			tot_len += 1
			for ind, geno in zip(inds, var[c][pos]):
				if geno == '-':
					indmiss[ind] += 1

	to_drop = []
	for ind in inds:
		if tot_len == 0:
			permiss = 1
		else:
			permiss = indmiss[ind] / float(tot_len)
		if permiss >= 0.8:
			to_drop.append(ind)
	
	return to_drop

def make_adegenet(inds, b, stem, drop):
	out = os.path.join(outdir, '%s.snp' % stem)
	o = open(out, 'w')
	o.write('>>>> begin comments - do not remove this line <<<<\n')
	o.write('>>>> end comments - do not remove this line <<<<\n')

	snp = {'0': {'0': '0', '1': '1'}, 
		   '1': {'0': '1', '1': '2'}, 
		   '-9': {'-9': '-'}}

	for ix, ind in enumerate(inds):
		if ind not in drop:
			gen = ''

			for c in sorted(b.keys()):
				for pos in sorted(b[c].keys()):
					geno = list(b[c][pos][ix])
					gen += snp[geno[0]][geno[1]]
			o.write('> %s \n%s\n' % (ind, gen))

	o.close()


def make_admixture(inds, var, stem, drop):
	out1 = os.path.join(outdir, '%s.ped' % stem)
	out2 = os.path.join(outdir, '%s.fam' % stem)
	o = open(out1, 'w')
	f = open(out2, 'w')

	for ix, ind in enumerate(inds):
		if ind not in drop:
			a = ['0', ind, '0', '0', '0', '0']
			f.write(' '.join(a) +' \n')
			for c in sorted(var.keys()):
				for pos in sorted(var[c].keys()):
					snp = code[var[c][pos][ix]]
					bp1 = snp[0]
					bp2 = snp[1]

					if bp1 == 'N':
						bp1 = '0'
					if bp2 == 'N':
						bp2 = '0'

					a.append(bp1)
					a.append(bp2)
			o.write(' '.join(a) + '\n')
	o.close()
	f.close()

	out3 = os.path.join(outdir, '%s.map' % stem)
	o = open(out3, 'w')
	tot_len = 0
	for c in var:
		pos = max([int(x) for x in var[c].keys()])
		tot_len += pos
	tot_len = (round(tot_len / 1000000, 2) * 100 + 1) * 10000

	cur_len = 0

	ix = 1
	for c in sorted(var.keys()):
		for pos in sorted(var[c].keys()):
			cur_len += int(pos)
			o.write('1\trs%s\t%.3f\t%s\n' % 
				(ix, cur_len / float(tot_len), cur_len))
			ix += 1

	o.close()


def get_sample(var, sample):
	sample = int(sample)
	var2 = {}
	while sample > 0:
		rc = random.choice(var.keys())
		if rc not in var2:
			var2[rc] = {}
		rcpos = random.choice(var[rc].keys())
		if rcpos not in var2[rc]:
			var2[rc][rcpos] = var[rc][rcpos]
			sample = sample - 1
	return var2 

MISS = float(args.MISS)
MAC = int(args.MAC)

code, invcode = make_code()

if args.drop:
	drop = re.split(',', args.drop)
else:
	drop = []

var, inds = parse_snps(snps, invcode, MISS, MAC, drop)
stem = args.stem + '.miss%s.MAC%s' % (MISS, MAC)

# randomly subsample one SNP per locus
if args.random:
	var = get_random(var)
	stem = stem + '.thinned'
if args.sample:
	var = get_sample(var, args.sample)
	stem = stem + '.sample%s' % args.sample

drop = remove_miss(var, inds)
if args.drop:
	drop = drop + re.split(',', args.drop)
drop = list(set(drop))
out = os.path.join(outdir, stem + '.drop.out')
of = open(out, 'w')
for l in drop:
	of.write(l + '\tdrop\n')

# get binary sructure
binary = make_binary(var, code)

if args.structure:
	make_structure(inds, binary, stem, code, drop)
if args.svd:
	make_svd(inds, var, stem, code, drop)
if args.adegenet:
	make_adegenet(inds, binary, stem, drop)
if args.admixture:
	make_admixture(inds, var, stem, drop)

of.close()
