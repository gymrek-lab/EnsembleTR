#!/usr/bin/env python3

import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import sys

scope = ['https://spreadsheets.google.com/feeds',
         'https://www.googleapis.com/auth/drive']

credentials = ServiceAccountCredentials.from_json_keyfile_name(
    'capillaryelectrophoresis-e683d5dd37a2.json', scope)
gc = gspread.authorize(credentials)
gsheet="1000GenomesRepeatValidationDatabase"

def LoadGSheet(wks):
    data = wks.get_all_values()
    headers = data.pop(0)
    return pd.DataFrame(data, columns=headers)

loci_data = LoadGSheet(gc.open(gsheet).worksheet("loci"))
primer_data = LoadGSheet(gc.open(gsheet).worksheet("primers"))
genotype_data = LoadGSheet(gc.open(gsheet).worksheet("genotypes"))
genotype_data = pd.merge(genotype_data, primer_data[["PrimerID","LocusID"]], on=["PrimerID"])
genotype_data = pd.merge(genotype_data, loci_data[["LocusID","Chrom","Start","Motif"]], on=["LocusID"])

for i in range(genotype_data.shape[0]):
    items = [genotype_data["Chrom"].values[i].strip(), \
             genotype_data["Start"].values[i].strip(), \
             genotype_data["LocusID"].values[i].strip(), \
             genotype_data["Motif"].values[i].strip(), \
             genotype_data["SampleID"].values[i].strip(), \
             genotype_data["Genotype"].values[i].strip()]
    if genotype_data["Genotype"].values[i].strip() == "": continue
    sys.stdout.write("\t".join([str(item) for item in items])+"\n")

