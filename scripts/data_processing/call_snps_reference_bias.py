import argparse
import os
import pandas as pd
import re
import subprocess
import gzip

"""
Sonal Singhal
created on 3 March 2023
Written assuming:
	* samtools 1.3.1
	* picard 2.4.1
	* ngm
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Reference correct genome. Written assuming "
				" samtools 1.3.1, picard 2.4.1, and"
				" bwa",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# sample
	parser.add_argument(
		'--lineage',
		type=str,
		default=None,
		help='Lineage for which to run script.'
		)

	# file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with data on individuals & species.'
		)

	# bwa
	parser.add_argument(
		'--bwa',
		type=str,
		default=None,
		help='bwa executable, full path.'
		)

	# samtools
	parser.add_argument(
		'--samtools',
		type=str,
		default=None,
		help='samtools executable, full path.'
		)

	# bcftools
	parser.add_argument(
		'--bcftools',
		type=str,
		default=None,
		help='bcftools executable, full path.'
		)

	# picard
	parser.add_argument(
		'--picard',
		type=str,
		default=None,
		help='picard executable, full path.'
		)
	
	# CPUs
	parser.add_argument(
		'--CPU',
		type=int,
		default=1,
		help='# of CPUs to use in alignment.'
		)

	# memory
	parser.add_argument(
		'--mem',
		type=int,
		default=1,
		help='Memory available, as an int, in terms of Gb.'
		)
		   
	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Output directory for alignments.'
		)

	# readdir
	parser.add_argument(
		'--readdir',
		type=str,
		default=None,
		help="Full path to files with reads."
		)

	# genome
	parser.add_argument(
		'--genome',
		type=str,
		default=None,
		help="Full path to  genome."
		)

	return parser.parse_args()


def get_info(args):
	# get the genome
	genome = args.genome

	outdir = args.outdir
	
	d = pd.read_csv(args.file)
	d = d.loc[d['ddRAD'] == True, ]
	d = d.loc[d['RECENT_RUN'] == False, ]
	samps = d.loc[d['CURRENT_TAXON'] == args.lineage, 'SAMPLE_ID'].tolist()
	reads = {}
	for s in samps:
		reads[s] = [os.path.join(args.readdir, '%s_R1_.fastq.gz' % s), 
				os.path.join(args.readdir, '%s_R2_.fastq.gz' % s)]

	return reads, genome, outdir


def prepare_seq(args, genome):
	# does all the prep necessary for the PRG
	if not os.path.isfile(genome + '.fai'):
		subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
	out = re.sub('.fa.*', '.dict', genome)
	if not os.path.isfile(out):
		subprocess.call("%s CreateSequenceDictionary R=%s O=%s" % 
				(args.picard, genome, out), shell=True)


def align_seq(args, sample, r, seq):
	dir = os.getcwd()

	out1 = os.path.join(dir, '%s.sam' % sample)
	out1s = os.path.join(dir, '%s.bam' % sample)
	out2 = os.path.join(dir, '%s.mateFixed.bam' % sample)
	out3 = os.path.join(dir, '%s.mateFixed.sorted.bam' %  sample)
	out4 = os.path.join(dir, '%s.rg.mateFixed.sorted.bam' % (sample))

	# need a tmpdir for when sorting BAM files
	tmpdir = os.path.join('/media/babs/brains3/', sample)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align using bwa because fastx
	out = seq + '.bwt'
	if not os.path.isfile(out):
		subprocess.call("%s index %s" % (args.bwa, seq), shell = True)
	subprocess.call("%s mem %s %s %s -o %s -t %s" % (args.bwa, seq, r[0], r[1], out1, args.CPU), shell=True)
	# fixmate
	# note that had used samtools to fixmates, appears to not work properly
	subprocess.call("%s view -b %s > %s" % (args.samtools, out1, out1s), shell=True)
	subprocess.call("%s fixmate %s %s" % (args.samtools, out1s, out2), shell=True)
	# sorted
	subprocess.call("%s sort -O bam -o %s --threads %s -m 2G -T %s %s" % (args.samtools, out3, args.CPU, tmpdir, out2), shell=True)
	# readgroup
	subprocess.call("%s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
			  (args.picard, out3, out4, sample, sample, sample), shell=True)
	# index "final" snp
	subprocess.call("%s index %s" % (args.samtools, out4), shell=True)
	
	# remove the files
	[os.remove(x) for x in [out1, out1s, out2, out3]]
	
	# remove the dir
	os.rmdir(tmpdir)

	return out4


def call_snps(args, files, genome):
	outdir = args.outdir

	out1 = os.path.join(outdir, '%s.bcf' % (args.lineage))
	out2 = os.path.join(outdir, '%s.vcf.gz' % (args.lineage))

	if not os.path.isfile(out2):
		subprocess.call("%s mpileup -A -B -I -O u -a FORMAT/DP --threads %s -f %s -o %s %s" % (args.bcftools, args.CPU, genome, out1, ' '.join(files)), shell=True)
		# output variant sites only
		subprocess.call("%s call -mO z -o %s %s" % (args.bcftools, out2, out1), shell=True)


def main():
	# get arguments
	args = get_args()
	reads, genome, outdir = get_info(args)

	inds = sorted(list(reads.keys()))
	# round 1
	# prep sequence
	prepare_seq(args, genome)
	# do the alignments
	bamfiles1 = []
	for ind in inds:
		bamout = align_seq(args, ind, reads[ind], genome)
		bamfiles1.append(bamout)
	genome1 = call_snps(args, bamfiles1, genome)
	xy = [os.remove(x) for x in bamfiles1]
	xy = [os.remove(x + '.bai') for x in bamfiles1]


if __name__ == "__main__":
	main()
