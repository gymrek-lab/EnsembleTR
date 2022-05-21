import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

pop = argv[1]
chr = argv[2]
samples = []

pop_dict = {"EUR":"European", "AFR": "African"}
pop_color = {'EUR':"blue", "AFR":"orange"}

with open(pop + "_samples.txt") as f:
    for line in f:
        samples.append(line.strip())


vcf_addr="/projects/ps-gymreklab/helia/ensembl/ensemble_out/merged_chr" +  chr  + "_sorted_ver2.vcf.gz"
all_imputed = pd.read_csv(vcf_addr, header=None,delim_whitespace=True, usecols=[1],comment='#')
all_gt = pd.read_csv(vcf_addr, header=None,delim_whitespace=True, usecols=[1],comment='#')
all_imputed.columns = ["pos"]
all_gt.columns = ["pos"]
all_imputed['pos'] = all_imputed['pos'].astype('int').astype('str')
all_gt['pos'] = all_gt['pos'].astype('int').astype('str')
all_gt = all_gt.drop_duplicates(subset='pos', keep="first")
all_imputed = all_imputed.drop_duplicates(subset='pos', keep="first")

for sample in samples:
    print(sample)
    df = pd.read_csv("diff_" + pop + "/" + chr + "/" + sample + "." + chr + ".diff.txt", header = None, delim_whitespace=True)
    df.columns = ['pos', 'gt1_' + sample, 'gt2_' + sample, 'i1_' + sample, 'i2_' + sample]
    df['gt_' + sample] = list(zip(df['gt1_' + sample], df['gt2_' + sample]))
    df['i_' + sample] = list(zip(df['i1_' + sample], df['i2_' + sample]))
    df['pos'] = df['pos'].astype('int').astype('str')
    df = df.drop_duplicates(subset='pos', keep="first")
    all_imputed = pd.merge(all_imputed, df[['pos', 'i_' + sample]], on='pos', how='left')
    all_gt = pd.merge(all_gt, df[['pos', 'gt_' + sample]], on='pos', how='left')

all_gt = all_gt.dropna(axis = 0, thresh=20)
all_imputed = all_imputed.dropna(axis = 0, thresh=20)

all_gt.to_csv("diff_" + pop + "/" + chr + "/" + "genotypes.csv", index=False)
all_imputed.to_csv("diff_" + pop + "/" + chr + "/" + "imputed.csv", index=False)

locus_dict = {}
for i in range(all_imputed.shape[0]):
    gts = list(all_gt.iloc[i])
    imputes = list(all_imputed.iloc[i])
    sum_11 = 0
    if gts[0] != imputes[0]:
        print("shit", gts[0])
    for j in range(1,len(samples)+1):
        sum_11 += int(gts[j] == imputes[j])
    locus_dict[gts[0]] = sum_11/len(samples)


plt.hist(locus_dict.values(), bins=60, color=pop_color[pop], edgecolor='#169acf', linewidth=0.5)
plt.xlabel("Concordance across " + str(len(samples)) + " samples of " + pop_dict[pop] + " population")
plt.ylabel("Number of loci")
plt.title("Chromosome " + chr)
plt.savefig("L1O_plots/" + "hist_" + pop + "_" + chr + ".jpg",dpi=1200)




