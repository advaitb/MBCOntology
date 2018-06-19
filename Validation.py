
# coding: utf-8

# Code for renaming processes in all the validation files

# In[20]:


import pandas as pd
import glob

class Validation:
    
    def customChange(self):
        #df = pd.read_csv('/Users/advaitbalaji/Desktop/name_change.txt', header = 0, sep = '\t')
        #old_p = df['Old Processes'].tolist()
        #new_p = df['New Processes'].tolist()
        pivot = 'transmembrane transport>'
        files = glob.glob('/Users/advaitbalaji/Desktop/Validation_results/*.txt')
        files.remove('/Users/advaitbalaji/Desktop/Validation_results/Validated_gene_scp_associations.txt')
        files.remove('/Users/advaitbalaji/Desktop/Validation_results/Validation_summary.txt')
        for file_r in files:
            w_f = '/Users/advaitbalaji/Desktop/NewValidation/'+file_r.split('/Users/advaitbalaji/Desktop/Validation_results/')[1]
            with open(file_r,'r') as f:
                new_lines = []
                lines = f.readlines()
                for line in lines:
                    if '@SCP:' in line:
                        #key  = line.split('@SCP: <')[1].split('>\t')[0]
                        #key  = line.split('@SCP: <')[1].split('\t')[0]
                        #if key in old_p:
                        if pivot in line:
                            #new_key = new_p[old_p.index(key)]
                            #new_key = key.split(pivot)[0]+'homeostasis>'
                            line = line.replace(pivot,'homeostasis>')
                            print('<------>DONE<------>')
                        new_lines.append(line)
                    else:
                        new_lines.append(line)
                with open(w_f,'w') as wf:
                    wf.writelines(new_lines)
            print('REWRITTEN FILE: '+file_r.split('/Users/advaitbalaji/Desktop/Validation_results/')[1]+ ' Index: '+ str(files.index(file_r))+'/'+str(len(files))+'...') 

    def rearrangeProcess(self):
        df = pd.read_csv('/Users/advaitbalaji/Desktop/name_change.txt', header = None, names = ['Old','New'] sep = '\t')
        data = pd.read_csv('/Users/advaitbalaji/Desktop/Supplementary Table S32 - gene-SCP associations_Gene_symbol.txt', header = 0, sep = '\t')
        data = data[data['ProcessLevel'] != -1][['ProcessName','Symbol']]
        old_p = df['Old'].tolist()
        new_p = df['New'].tolist()
        files = glob.glob('/Users/advaitbalaji/Desktop/Validation_results/*.txt')
        genes_intersect = {}
        print('Creating Dict!...')
        for p in old_p:
            old_g = data[data['ProcessName'] == p]['Symbol'].unique().tolist()
            new_g = data[data['ProcessName'] == new_p[old_p.index(p)]]['Symbol'].unique().tolist()
            i_g = [g for g in old_g if g in new_g]
            genes_intersect[p] = i_g
        flag = False
        #files.remove('/Users/advaitbalaji/Desktop/Validation_results/Validated_gene_scp_associations.txt')
        #files.remove('/Users/advaitbalaji/Desktop/Validation_results/Validation_summary.txt')
        for file_r in files:
                w_f = '/Users/advaitbalaji/Desktop/NewValidation/'+file_r.split('/Users/advaitbalaji/Desktop/Validation_results/')[1]
                with open(file_r,'r') as f:
                    new_lines = []
                    lines = f.readlines()
                    for line in lines:
                        if '---------------------------------END------------------------' or '<--------------------------------------------------------------------------------------------------------------------------------------->' in line:
                            flag = False
                            new_lines.append(line)
                        elif '@SCP:' in line:
                            new_lines.append(line)
                            key  = line.split('@SCP: <')[1].split('>\t')[0]
                            gene = line.split('\t: <')[1].split('> -')[0]
                            if key in old_p:
                                genes = genes_intersect[key]
                                if genes == []:
                                    continue
                                if gene in genes:
                                    flag = True
                        elif '+Validation: <T>' in line and flag:
                            line.replace('+Validation: <T>','+Validation: <C>')
                            new_lines.append(line)
                        else:
                            new_lines.append(line)
                    with open(w_f,'w') as wf:
                        wf.writelines(new_lines)
                print('REWRITTEN FILE: '+file_r.split('/Users/advaitbalaji/Desktop/Validation_results/')[1]+ ' Index: '+ str(files.index(file_r))+'/'+str(len(files))+'...')

if __name__ == '__main__': 
    
    v = Validation()
    v.rearrangeProcess()         
            
            
            
            
            