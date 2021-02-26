#!/usr/bin/env python3

"""
Usage:

./convert_cap.py <capfile> <refvcf>

Outputs capillary data in a VCF format that will be comparable to refvcf
"""

import sys
import vcf
import pandas as pd

try:
    capfile = sys.argv[1]
    callvcf = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

# Get sample list
capdata = pd.read_csv(capfile, sep="\t", names=["chrom","pos","LocusID","Motif","sample","genotype"])
samples = list(set(capdata["sample"]))

# Write VCF header
sys.stdout.write('##fileformat=VCFv4.1\n')
sys.stdout.write('##command=GangSTR dummy\n')
sys.stdout.write('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">\n')
sys.stdout.write('##INFO=<ID=CRU,Number=1,Type=String,Description="Capillary repeat motif">\n')
sys.stdout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
sys.stdout.write('##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period (length of motif)">\n')
sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
sys.stdout.write('##FORMAT=<ID=REPCN,Number=2,Type=Integer,Description="Genotype given in number of copies of the repeat motif">\n')
sys.stdout.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+samples)+"\n")

def GetAlt(count, motif):
    return motif.upper()*count

def GetAlleleDict(refgt, capgts, motif):
    # Get list of alts
    alts = []
    for i in range(len(capgts)):
        try:
            a1, a2 = [int(item) for item in capgts[i].strip().split(",")]
        except ValueError: continue
        if a1 != refgt: alts.append(a1)
        if a2 != refgt: alts.append(a2)
    alts = sorted(list(set(alts)))
    
    # Map capillary gt -> alt num
    adict = {refgt: 0}
    for i in range(len(alts)):
        adict[alts[i]] = i+1
    return [GetAlt(item, motif) for item in alts], adict

# Match each locus to reference VCF
vcf_reader = vcf.Reader(open(callvcf, "rb"))
W = 10 # search +/- this many bp
caploci = capdata[["chrom","pos","Motif"]].drop_duplicates()
for i in range(len(caploci)):
    cap_chrom = caploci["chrom"].values[i]
    cap_pos = caploci["pos"].values[i]
    cap_motif = caploci["Motif"].values[i]
    try:
        records = vcf_reader.fetch(cap_chrom, cap_pos-W, cap_pos+W)
    except ValueError:
        sys.stderr.write("No match for %s:%s:%s\n"%(cap_chrom, cap_pos, cap_motif))
        continue
    for r in records:
        #sys.stderr.write("Found match for %s:%s:%s\n"%(cap_chrom, cap_pos, cap_motif))
        #sys.stderr.write(str(r) + "\n")
        if len(r.INFO["RU"])!=len(cap_motif):
            sys.stderr.write("Motifs don't match. skipping %s:%s:%s\n"%(cap_chrom, cap_pos, cap_motif))
            continue
        # Set up info
        info_strings = ["RU=%s"%r.INFO["RU"], \
                        "CRU=%s"%cap_motif, \
                        "END=%s"%r.INFO["END"], \
                        "PERIOD=%s"%r.INFO["PERIOD"]]
        
        # Set up alts
        refgt = int(len(r.REF)/r.INFO["PERIOD"])
        alts, adict = GetAlleleDict(refgt, list(capdata[(capdata["chrom"]==cap_chrom) & (capdata["pos"]==cap_pos)]["genotype"]), r.INFO["RU"])

        # Get sample info
        sample_info = []
        for s in samples:
            x = capdata[(capdata["chrom"]==cap_chrom) & (capdata["pos"]==cap_pos) & (capdata["sample"]==s)]
            if x.shape[0] == 0:
                sample_info.append(".")
                continue
            if x.shape[0] > 1:
                sys.stderr.write("Error found multiple capillary entries for %s for the same locus\n"%s)
                sys.exit(1)
            try:
                a1,a2 = [int(item) for item in x["genotype"].values[0].strip().split(",")]
            except:
                sys.stderr.write("Could not parse genotype %s\n"%x["genotype"].values[0].strip())
                continue
            sample_gt = "%s/%s"%(adict[a1],adict[a2])
            sample_repcn= "%s,%s"%(a1, a2)
            sample_info.append("%s:%s"%(sample_gt, sample_repcn))
        # Output record should match chrom,pos,repeat unit
        altfield = "."
        if len(alts)>0: altfield = ",".join(alts)
        items = [cap_chrom, r.POS, ".", r.REF, altfield, ".", ".", ";".join(info_strings), "GT:REPCN"]+sample_info
        sys.stdout.write("\t".join([str(item) for item in items])+"\n")


