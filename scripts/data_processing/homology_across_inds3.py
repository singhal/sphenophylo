import argparse
import subprocess
import glob
import os
import multiprocessing as mp

dir = '/home/babs/spheno/phylogeny/ddRAD/alignments/'

parser = argparse.ArgumentParser(description='align a directory.')
parser.add_argument('--n', help="cluster name")
args = parser.parse_args()

subdir = os.path.join(dir, args.n)
alns = glob.glob(subdir + '/*fa')

def run_aln(alnfile):
    subprocess.call("mafft --adjustdirection %s > %s.aln" % (alnfile, alnfile), shell = True)

pool = mp.Pool(24)
alnruns = pool.map(run_aln, alns)

