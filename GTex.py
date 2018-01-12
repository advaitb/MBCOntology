import pandas as pd
import csv

#Extract symbols for background

df = pd.read_csv('/Users/advaitbalaji/Desktop/gene_all_info_Homo_sapiens_2016June23_2016June23.jns', sep = '\t', header = 0)
df = df[df.EnsemblGene != 'E_m_p_t_y']
df.reset_index(inplace = True)
del df['index']
#print(df.head())
eg = df.EnsemblGene.tolist()
#print(eg)
sym = df.Symbol.tolist()
#print(sym)
#print('\n')
#print(len(eg)),print(len(sym))
#if len(eg) == len(sym):
 #   print('Yippee!')
#else:
 #   print('Bleh!')


with open('/Users/advaitbalaji/Desktop/background_genes.txt','r') as inp, open('/Users/advaitbalaji/Desktop/background_genes_correct.txt','w') as out:
    lines = inp.readlines()
    for line in lines:
        line = line.split('\n')[0]
        if line in eg: 
            temp_sym = sym[eg.index(line)].upper()
            print("Writing symbol "+ temp_sym+"!")
            out.write(temp_sym+'\n')

  

#Write new GCT file with symbols

with open('/Users/advaitbalaji/Desktop/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct', 'r') as initf, open('/Users/advaitbalaji/Desktop/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm_FINALAB.gct','a') as output:
    initlines = initf.readlines()
    header, lines = initlines[0].split('\t'), initlines[1:]
    #writer = csv.writer(output, delimiter = '\t')
    #writer.writerow(header)
    output.write('\t'.join(header))
    for line in lines:
        #print(line)
        #print('\n')
        #print(line.split('\t')[0])
        temp =  line.split('\t')
        ens_gene = temp[0].split('.')[0]
        #print(ens_gene)        
        if ens_gene in eg:
            sym_gene = sym[eg.index(ens_gene)]
            #print(sym_gene)
            #print('\n')
            temp[0] = sym_gene.upper() 
            #print(temp)  
            output.write('\t'.join(temp))


#write CSV with top 500 hits for each GTEx Tissue

import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns

import numpy as np
from scipy import stats
import math
import csv
import pandas as pd
np.seterr(divide='ignore', invalid='ignore')
with open('/Users/advaitbalaji/Desktop/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm_FINAL.gct', 'r') as f:
    lines = f.readlines()
    keys, genes = lines[0].split('\t')[2:], lines[1:]
    z_score = {key: [] for key in keys}
    #print(header)
    for line in genes:
        #gene_list.append(line.split('\t')[0].split('.')[0])
        #print(gene_list)
        #break
        sentence = stats.zscore(np.array([float(i) for i in line.split('\t')[2:]]))
        if math.isnan(max(sentence)):
            continue
        p = [(line.split('\t')[0].split('.')[0], value) for value in np.array(sentence).tolist()]
        for (key,tup) in zip(z_score,p):
            z_score[key].append(tup)
    #print(z_score)
    #print('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
    for key in z_score:
        z_score[key].sort(key = lambda x: x[1], reverse = True)
        z_score[key] = z_score[key][:500]
        for i in range(len(z_score[key])):
            a,b = z_score[key][i]
            z_score[key][i] = a
    keys = sorted(z_score.keys())
    with open("/Users/advaitbalaji/Desktop/output.csv", "a") as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(z_score.keys())
        for row in zip(*z_score.values()):
            writer.writerow(list(row))
    #print(z_score)

            
#Write background genes            

    with open('/Users/advaitbalaji/Desktop/background_genes.txt', 'a') as of:
        for g in gene_list:
            of.write(g+'\n')
'''
with open('/Users/advaitbalaji/Desktop/background_genes.txt','r') as inf, open('/Users/advaitbalaji/Desktop/background_genes_f.txt','a') as outf:
    lines = inf.readlines()
    for line in lines:
        outf.write(line.split('.')[0]+'\n')
'''