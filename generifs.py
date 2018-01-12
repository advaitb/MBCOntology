import pandas as pd
import numpy as np
from statistics import mean,median
import time
import dask.dataframe as dd
import sys
import sqlite3
from collections import Counter
import csv
'''
df = pd.read_csv('/Users/advaitbalaji/Desktop/Generifs9606.csv', sep = '\t', header = 0)
genes = df['Gene ID'].unique().tolist()
abstract_count = []
for gene in genes:
    count = len(df[df['Gene ID'] == gene]['PubMed ID (PMID) list'].unique().tolist())
    abstract_count.append(count)
print("Maximum number of abstracts/gene: "+ str(max(abstract_count))+'\n')
print("Minimum number of abstracts/gene: "+ str(min(abstract_count))+'\n')
print("Median number of abstracts/gene: "+ str(median(abstract_count))+'\n')
print("Average number of abstracts/gene: "+ str(mean(abstract_count))+'\n')
print(len(abstract_count))
'''

'''
df = pd.read_csv('/Users/advaitbalaji/Downloads/AdvaitB/Suppl. table 32 - mbc_association_pvalue_rank_Gene_symbol_sameLevel_pvalue_sameChildrenSet_pvalue.txt', sep ='\t', header = 0)
back_df = df[df.ProcessLevel == -1]
total_genes = back_df.Symbol.unique().tolist()
data = pd.read_csv('/Users/advaitbalaji/Desktop/Test_gtex/gene_all_info_Homo_sapiens_2016June23_2016June23.jns', sep = '\t', header = 0)
sym = data.Symbol.unique().tolist()
gid = data.GeneID.unique().tolist()

with open('/Users/advaitbalaji/Desktop/Genesymback.txt','w') as f:
    f.write('Symbol,GeneID\n')
    
    for gene in total_genes:
        if gene in sym:
            index = sym.index(gene)
            gene_id = gid[index]
            f.write(gene+','+str(gene_id)+'\n')
'''
'''
data = pd.read_csv('/Users/advaitbalaji/Desktop/Generifs9606.csv', sep = '\t', header = 0)
#df = pd.read_csv('/Users/advaitbalaji/Desktop/Genesymback.txt', sep = ',', header = 0)
#geneids = df.GeneID.unique().tolist()
#arr = []
#for gid in geneids:
#    arr.append(data[data['Gene ID'] == gid])
#final_frame = pd.concat(arr, ignore_index = True)

#print(final_frame.head(20))
#print("Finished background genes containing final frame!")
#(final_frame.head(20))
#print(len(final_frame['Gene ID'].unique().tolist()))
#abstracts = final_frame['PubMed ID (PMID) list'].unique().tolist()
abstracts = data['PubMed ID (PMID) list'].unique().tolist()
print(len(abstracts))
print("Abstracts Ready...")
#print(len(abstracts))
#ab_dict = {}

with open('/Users/advaitbalaji/Desktop/ab_pairs.txt','a+') as f:
    f.write('Abstract1'+'\t'+'Abstract2'+'\t'+'Weights'+'\n')
    for i in range(0,46):
        abstract = abstracts[10000*i:10000+10000*i]
        modified_content = [(abstract[k],abstract[l]) for k in range(len(abstract)) for l in range(k+1, len(abstract))]
        for pair in modified_content:
            a,b = pair
            f.write(str(a)+'\t'+str(b)+'\t'+str(0)+'\n')
        print("Dictionary "+str(i+1)+" done...")

with open('/Users/advaitbalaji/Desktop/ab_pairs.txt','a+') as f:
    abstract = abstracts[450000:451847]
    modified_content = [(abstract[k],abstract[l]) for k in range(len(abstract)) for l in range(k+1, len(abstract))]
    for pair in modified_content:
        a,b = pair
        f.write(str(a)+'\t'+str(b)+'\t'+str(0)+'\n')
    print("Dictionary 47 done...")
    print('File Written Successfully!')
'''
'''
start_time = time.time()
with open('/Users/advaitbalaji/Desktop/test.txt','r') as f:
    lines = f.readlines()
    print(lines[1])
    print(lines[500005])
'''
#pd.to_pickle('/Users/advaitbalaji/abstract_data.pkl')
#print(time.time()-start_time)
#sys.exit()
#print('\n')
#print(df.head())



connex = sqlite3.connect("abstract.db")  
cur = connex.cursor()
print('CONNECTED!...')
sql = "SELECT * FROM abstract WHERE Weights != 0"
head = True
'''
for chunk in pd.read_csv("/Users/advaitbalaji/Desktop/ab_pairs.txt", chunksize=1000000, sep = '\t'):
    chunk.to_sql(name="abstract", con=connex, if_exists="append", index=False)
print("CREATED DATABASE SUCCESFULLY!...")   
cur.execute('CREATE INDEX indx ON abstract(Abstract1, Abstract2)')
print("CREATED INDEX SUCCESSFULLY!...")
connex.commit()
connex.close()
'''

'''
df = pd.read_csv('/Users/advaitbalaji/Desktop/Generifs9606.csv', sep = '\t', header = 0)
genes = df['Gene ID'].unique().tolist()
for gene in genes:
    abstracts = df['PubMed ID (PMID) list'][df['Gene ID'] == gene].unique().tolist()
    modified_content = [(abstracts[k],abstracts[l]) for k in range(len(abstracts)) for l in range(k+1, len(abstracts))]
    print("PAIRS CREATED!...")
    for pair in modified_content:
        a,b = pair
        #print(a)
        #print(b)
        #cur.execute("SELECT Weights FROM abstract INDEXED BY indx WHERE Abstract1 = %s AND Abstract2 = %s" %(str(a),str(b)))
        #p = cur.fetchone()
        #a =p[0] + 1
        #print(p[0]+1)
        #print('\n')
        cur.execute("UPDATE abstract SET Weights = Weights + 1 WHERE Abstract1 = %s AND Abstract2 = %s" %(str(a),str(b)))
        connex.commit()
        #p = cur.fetchone()
        #print(p)
        #connex.commit()
        #print("IT WORKED!")
    print("Gene "+ str(genes.index(gene))+"/"+str(len(genes))+" completed...")
'''
'''
df = pd.read_csv('/Users/advaitbalaji/Desktop/p.txt', sep = '\t', header = 0)
df.head()
connex = sqlite3.connect('database.db')
cur = connex.cursor()
df.to_sql(name = 'database', con = connex, if_exists = 'replace', index = False)
'''
with open('/Users/advaitbalaji/Desktop/Graph.csv','a+') as f:
    for chunk in pd.read_sql(sql, connex, chunksize = 1000000):
        if head:
            chunk.to_csv(f, sep = '\t', header = head , index = False)
            print("WRITTEN CHUNK SUCCESSFULLY!...")
            head = False        
        else:
            chunk.to_csv(f, sep = '\t', header = head, index = False)
            print("WRITTEN CHUNK SUCCESSFULLY!...")
    print("WRITTEN CSV SUCCESSFULLY!...")
            
'''  
with open('/Users/advaitbalaji/Desktop/graph.csv','w') as f:
    writer = csv.writer(f, delimiter = '\t')
    writer.writerow(['Abstract1','Abstract2','Weights'])
    cur.execute("SELECT * FROM abstract")
    k = cur.fetchall()
    for trip in k:
        a,b,c = trip
        writer.writerow([a,b,c])
        print('WRITING INTO CSV SUCCESSFUL!...')
'''    

    
'''
    with open('/Users/advaitbalaji/Desktop/test.txt','a+') as f:
        f.write('Abstract1'+'\t'+'Abstract2'+'\t'+'Weights')
        for key in ab_dict:
            a,b = key
            f.write(str(a)+'\t'+str(b)+'\t'+str(0)+'\n')
'''
'''
modified_content = [(abstracts[k],abstracts[l]) for k in range(len(abstracts)) for l in range(k+1, len(abstracts))]
print("Finished creating abstract gene pairs!....")
abstract_count = {key: 0 for key in modified_content}
print("Abstract key successfully created! PHEW!")
'''

'''
for i in range(0,33):
    
    abstracts = abstracts[10000*i:10000+10000*i]


    modified_content = [(abstracts[k],abstracts[l]) for k in range(len(abstracts)) for l in range(k+1, len(abstracts))]
    #with open('/Users/advaitbalaji/Desktop/GenepairsNov2017.txt', 'w') as wf:
    #    for i in range(len(modified_content)):
    #        a,b = modified_content[i]
    #        wf.write(a+','+b+'\n')
    #print('Created gene pairs...')
    #def commonelements( a, b):
    #    return list(set(a) & set(b))
    print("Finished Creating abstract pairs! Phew!")
    print("Starting Abstract Count....")
    with open('/Users/advaitbalaji/Desktop/CommonGenes.txt', 'a+') as f:   
        f.write('Source,Target, Weights'+'\n')
        for pair in range(len(modified_content)):
            abstracta,abstractb = modified_content[pair]
            genesa = final_frame[final_frame['PubMed ID (PMID) list'] == abstracta]['Gene ID'].unique().tolist()
            genesb = final_frame[final_frame['PubMed ID (PMID) list'] == abstractb]['Gene ID'].unique().tolist()
            #totala = len(pgenesa)
            #totalb = len(pubmedb)
            common_genes = [gene for gene in genesa if gene in genesb]
            common_count = len(common_genes)
            #print(pp)
            print('Scanned abstract pair '+ str(pair+1)+'/'+str(len(modified_content)))
            #if pp <= 1:
            #    continue
            #pnp = totalb - pp
            #npp = totala - pp
            #temp = 662917 - totala
            #npnp = temp-pnp
            #oddsratio, pvalue =  stats.fisher_exact([[pp,pnp],[npp,npnp]])
            if common_count == 0:
                continue
            f.write(str(abstracta)+','+str(abstractb)+','+str(common_count)+'\n')


abstracts = abstracts[330000:337769]


modified_content = [(abstracts[i],abstracts[j]) for i in range(len(abstracts)) for j in range(i+1, len(abstracts))]
#with open('/Users/advaitbalaji/Desktop/GenepairsNov2017.txt', 'w') as wf:
#    for i in range(len(modified_content)):
#        a,b = modified_content[i]
#        wf.write(a+','+b+'\n')
#print('Created gene pairs...')
#def commonelements( a, b):
#    return list(set(a) & set(b))
print("Finished Creating abstract pairs! Phew!")
print("Starting Abstract Count....")
with open('/Users/advaitbalaji/Desktop/CommonGenes.txt', 'a') as f:   
    f.write('Source,Target, Weights'+'\n')
    for i in range(len(modified_content)):
        abstracta,abstractb = modified_content[i]
        genesa = final_frame[final_frame['PubMed ID (PMID) list'] == abstracta]['Gene ID'].unique().tolist()
        genesb = final_frame[final_frame['PubMed ID (PMID) list'] == abstractb]['Gene ID'].unique().tolist()
        #totala = len(pgenesa)
        #totalb = len(pubmedb)
        common_genes = [gene for gene in genesa if gene in genesb]
        common_count = len(common_genes)
        #print(pp)
        print('Scanned abstract pair '+ str(i+1)+'/'+str(len(modified_content)))
        #if pp <= 1:
        #    continue
        #pnp = totalb - pp
        #npp = totala - pp
        #temp = 662917 - totala
        #npnp = temp-pnp
        #oddsratio, pvalue =  stats.fisher_exact([[pp,pnp],[npp,npnp]])
        f.write(str(abstracta)+','+str(abstractb)+','+str(common_count)+'\n')
'''    