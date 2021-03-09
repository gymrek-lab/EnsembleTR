#!/usr/bin/env python3

"""
Usage:

./combine_freqs.py <freqfile>

Combines frequencies from each population group to a single output file for plotting
"""

import sys
import pandas as pd

try:
    freqfile = sys.argv[1]
    header = sys.argv[2].split(',')
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# get labels/sample size for each pop
POPS = [   "ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", \
           "CLM", "MXL", "PEL", "PUR", \
           "CDX", "CHB", "CHS", "JPT", "KHV", \
           "CEU", "FIN", "GBR", "IBS", "TSI", \
           "BEB", "GIH", "ITU", "PJL", "STU"]
ssize = {}
for pop in POPS: ssize[pop] = 100 # TODO!!! change to actual sample sizes...

popinfo = {}
for pop in ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"]: popinfo[pop] = "AFR"
for pop in ["CLM", "MXL", "PEL", "PUR"]: popinfo[pop] = "AMR"
for pop in ["CDX", "CHB", "CHS", "JPT", "KHV"]: popinfo[pop] = "EAS"
for pop in ["CEU", "FIN", "GBR", "IBS", "TSI"]: popinfo[pop] = "EUR"
for pop in ["BEB", "GIH", "ITU", "PJL", "STU"]: popinfo[pop] = "SAS"

def GetCounts(x):
    n = ssize[x["pop"]]*2
    freqs = x["afreq"]
    counts = []
    for allele in freqs.split(","):
        al, fq = allele.split(":")
        ct = int(float(fq)*n)
        counts.append("%s:%s"%(al, ct))
    return ",".join(counts)

def AggCounts(countlist):
    countdict = {} # allele -> count
    for cl in countlist:
        for allele in cl.split(","):
            a1, ct = allele.split(":")
            countdict[a1] = countdict.get(a1, 0) + int(ct)
    total = float(sum(countdict.values()))
    freqs = []
    for item in countdict.keys():
        freqs.append("%s:%.3f"%(item, countdict[item]/total))
    return ",".join(freqs)

# Load data
# changed old names=["pop","chrom","start","end","freqs","het","numalleles"]
# with a more versatile header that updates based on content of stats file
data = pd.read_csv(freqfile, sep="\t", names=["pop"] + header)
data["superpop"] = data["pop"].apply(lambda x: popinfo[x])
data["counts"] = data.apply(GetCounts, 1)

# Process each locus
superpops = ["EUR","EAS","AFR","AMR","SAS"]
sys.stdout.write("\t".join(["chrom","pos"]+["freq_%s"%sp for sp in superpops])+"\n")
loci = data[["chrom","start"]].drop_duplicates()
for i in range(loci.shape[0]):
    chrom = loci["chrom"].values[i]
    pos = loci["start"].values[i]
    aggdata = data[(data["chrom"]==chrom) & (data["start"]==pos)].groupby(["superpop"], as_index=False).agg({"counts": AggCounts})
    freqlist = []
    for sp in superpops:
        freqlist.append(aggdata[aggdata["superpop"]==sp]["counts"].values[0])
    sys.stdout.write("\t".join([chrom, str(pos)]+freqlist)+"\n")
    



