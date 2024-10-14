import re
import os
import random
import subprocess
import glob
import argparse

parser = argparse.ArgumentParser(description="create pop gen runs")
parser.add_argument('--bed',
			default = None,
			help = "bedfile")
args = parser.parse_args()

# do admixture
file = re.sub('\.bed', '.ped', args.bed)
ind = 0
o = open(file, 'r')
for l in o:
	ind += 1
o.close()

base = re.sub(".*\/", "", args.bed)
base = re.sub("\.bed", "", base)
out = re.sub(".miss.*", "", base)

if ind > 2:
	maxk = ind
	if maxk > 13:
		maxk = 13
	for k in range(1, maxk):
                print(base, k)
                output = "/media/brains1/spheno/fastStructure/%s.%s.meanQ" % (out, k)
                if not os.path.isfile(output):
                        subprocess.call("~/bin/fastStructure -K %s --input=%s --output=%s --cv=10" % (k, base, out), shell = True)
