#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 12:31:49 2017

@author: advaitbalaji
"""
import pandas as pd
import pickle
import time
from itertools import chain
import numpy as np
'''

import pandas as pd
import networkx as nx
import numpy as np

df =pd.read_csv('/Users/advaitbalaji/Desktop/Graph.csv', sep = '\t', header =0)
print("LOADED DATAFRAME...")
G = nx.from_pandas_dataframe(df, source = 'Abstract1', target = 'Abstract2', edge_attr = 'Weights')
print("CREATED GRAPH...")
nx.write_gpickle(G, '/Users/advaitbalaji/Desktop/pickled_ab')
print("PICKLED AND DONE!...")
'''
'''
df = pd.read_csv('/Users/advaitbalaji/Downloads/MBC_ontology/Pubmed_search/Xml/MBC_ontology/Suppl. table 32 - mbc_association_pvalue_rank_Gene_symbol_sameLevel_pvalue_sameChildrenSet_pvalue.txt', sep = '\t', header = 0)
df = df[df.ProcessLevel == -1]
mbc_genes = df.Symbol.unique().tolist()


data = pd.read_csv('/Users/advaitbalaji/Downloads/generifs_basic', sep = '\t', header = 0)
data = data[data['#Tax ID'] == 9606]
#geneid = data['Gene ID'].unique().tolist()
data = data[['Gene ID','PubMed ID (PMID) list']]

#generif_genes = data['Gene ID'].unique().tolist()

#absent_genes = [gene for gene in generif_genes if gene not in mbc_genes]

#print(len(absent_genes))
#time.sleep(30)
temp = pd.read_csv('/Users/advaitbalaji/Downloads/Test_gtex/gene_all_info_Homo_sapiens_2016June23_2016June23.jns', sep = '\t', header =0)
temp = temp[['Symbol','GeneID']]

map_dict = {}

for gid in temp.GeneID.unique().tolist():
    symbol = ''.join(temp[temp.GeneID == gid].Symbol.unique().tolist())
    map_dict[gid] = symbol
    print("Gene "+str(temp.GeneID.unique().tolist().index(gid) + 1) +"/"+str(len(temp.GeneID.unique().tolist())) + " done...")


print(len(map_dict))
time.sleep(20)

#common_gid = [gid for gid in list(map_dict.keys()) if gid in data['Gene ID'].unique().tolist()]
keys = list(map_dict.keys())

final_dict = {'Abstract' : [],'GeneSymbol': []}


#with open("/Users/advaitbalaji/Desktop/missinggeneabstract.csv",'w') as f:
for pmid in data['PubMed ID (PMID) list'].unique().tolist():
    genes_per_abstract = []
    for geneid in data[data['PubMed ID (PMID) list'] == pmid]['Gene ID'].unique().tolist():
        if geneid in keys:
            if map_dict[geneid] not in mbc_genes:
                genes_per_abstract.append(map_dict[geneid])
    if genes_per_abstract:
        final_dict['Abstract'].append(pmid)
        final_dict['GeneSymbol'].append(genes_per_abstract)
    print("WRITTEN ABSTRACT "+str(data['PubMed ID (PMID) list'].unique().tolist().index(pmid) + 1) +"/"+str(len(data['PubMed ID (PMID) list'].unique().tolist()))+"...")
print(len(final_dict))
with open("/Users/advaitbalaji/Desktop/missinggeneabstract.pickle","wb") as f:
    pickle.dump(final_dict,f,protocol=pickle.HIGHEST_PROTOCOL)
    
#with open("/Users/advaitbalaji/Desktop/Gene-ID.txt","wb") as of:
#    pickle.dump(map_dict, of)
'''
'''
#import pickle
with open('/Users/advaitbalaji/Desktop/missinggeneabstract.pickle','rb') as rf:
    final_dict = pickle.load(rf)
print("MAKING PAIRS!...")
pmids = final_dict['Abstract']
genes_list = final_dict['GeneSymbol']
pmid_pairs = [(pmids[i],pmids[j]) for i in range(len(pmids)) for j in range(i+1,len(pmids))]
count = 0
print("MADE ALL PAIRS!...")
print("\n")
with open('/Users/advaitbalaji/Desktop/initial_abstract_check_rerun.csv','w') as of:
    of.write("Source,Target,Weights\n")
    for i in range(len(pmid_pairs)):
        a,b = pmid_pairs[i]
        genes_a = genes_list[pmids.index(a)]
        genes_b = genes_list[pmids.index(b)]
        weight = len([gene for gene in genes_a if gene in genes_b])
        if(weight >= 1):
            of.write(str(a)+","+str(b)+","+str(weight)+'\n')
            count+=1
        print("SCANNED GENE: "+str(i+1)+"/"+str(len(pmid_pairs))+"...")
'''

with open('/Users/advaitbalaji/Desktop/geneandnames.csv','w') as wf:
    wf.write('Symbol\tName\tSynonyms\n')
    df = pd.read_csv("/Users/advaitbalaji/Desktop/gene_list.txt", sep = '\t', header = 0)
    data = pd.read_csv("/Users/advaitbalaji/Desktop/geneids.csv", header = 0)
    S_list = df['Approved Symbol'].unique().tolist()
    for gene in data.GeneSymbol.unique().tolist():
        if gene in S_list:
            name = ''.join(df[df['Approved Symbol'] == gene]['Approved Name'].tolist())
            if np.nan in df[df['Approved Symbol'] == gene]['Synonyms'].tolist():
                synonym = ''
            else:
                synonym = ''.join(df[df['Approved Symbol'] == gene]['Synonyms'].tolist())
            wf.write(gene+'\t'+name+'\t'+str(synonym)+'\n')
        else:
            wf.write(gene+'\t '+'\t \n')
        print("WRITING!...")
