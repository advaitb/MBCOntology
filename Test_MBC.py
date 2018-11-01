import pandas as pd
from collections import OrderedDict
import pickle
from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Test:
    
    path = '/Users/advaitbalaji/'
    downloads = 'Downloads/'
    desktop = 'Desktop/'
    
    
    def read(self,destination,fname,delimiter,header):
        if destination.lower() == 'downloads':
            return pd.read_csv(self.path+self.downloads+fname, sep = delimiter , header = header)
        return pd.read_csv(self.path+self.desktop+fname, sep = delimiter , header = header)
    
    def openFileLineWise(self,fname):
        with open(self.path+self.downloads+fname,"r") as f:
            return f.readlines()
    
    def writeDictToFileLineWise(self, d, name):
        with open(self.path+self.desktop+name+'.txt','w') as outf:
            for i, j in d.items():
                outf.write(i+'\t'+' '.join(j)+'\n')
        print("WRITTEN!...")    
    
    
    
    def prepareOrderedDict(self,information):
        return{info.split('\t')[0].split('(GO')[0].strip(): [ a.split(',')[0] for a in info.split('\t')[2:-1] ] for info in information}
        
                 
    
    def orderDictionaryWithKey(self,d):
        return OrderedDict(sorted(d.items(), key = lambda t: t[0].lower())) 


    def pickleDictToBinary(self,destination,name,d): 
        print("PICKLED!...")
        if destination.lower() == 'desktop':
            with open(self.path+self.desktop+name+'.pickle',"wb") as wf:
                pickle.dump(d, wf, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(self.path+self.downloads+name+'.pickle',"wb") as wf:
                pickle.dump(d, wf, protocol=pickle.HIGHEST_PROTOCOL)
        
            
    def readPickledDictFromBinary(self,destination,name):
        print("UNPICKLED!...")
        if destination.lower() == 'desktop':
            with open(self.path+self.desktop+name+'.pickle',"rb") as rf:
                return pickle.load(rf)
        else:
            with open(self.path+self.downloads+name+'.pickle',"rb") as rf:
                return pickle.load(rf)
        

    def similar(self, a, b):
        print(SequenceMatcher(None, a, b).ratio())
        
    def optionFoRIU(self,option):
        if option.lower() == 'every':
            return False,True
        elif option.lower() == 'atleastone':
            return True,False
        else:
            return True,True
    
    def createMappingFromBackbone(self,directory,inputfile):
        with open(self.path+directory+inputfile,'r') as inpf,open(self.path+self.desktop+"Mapping_preliminary.txt",'w') as of:
            of.write("MBCO\tGO\n")
            lines = inpf.readlines()
            for line in lines:
                line = line.strip('-')
                try:
                    a,b = line.split('\t')[0],line.split('\t')[1]
                except IndexError:
                    continue
                of.write(a+'\t'+b)    
            print("PRELIMINARY MAPPING CREATED SUCCESSFULLY! CHECK FOR 3 COLUMNS, MAY NOT BE SEPARATED!...")
            print("USE CLEAN NaN TO AVOID Nan ERROR AND GENERATE FINAL MAPPING!...")
            
    def createFinalMappingWithoutNaN(self, directory, inputfile):
        df = pd.read_csv(self.path+directory+inputfile, sep = '\t', header = 0)
        df.dropna(axis = 0,inplace  = True)
        df = df.reset_index(drop=True)
        df.to_csv(self.path+self.desktop+'Mapping.txt', sep = '\t', header = True, index = False)
    
    def prepareGODictFromJensHansenFormat(self,directory,inputfile):
        df = pd.read_csv(self.path+directory+inputfile, sep = '\t', header = 0)
        df = df[['Source_name','Target']]
        process_names = df['Source_name'].unique().tolist()
        with open(self.path+directory+"GO_dict.txt","w") as of:
            of.write("Process\tGenes\n")
            for process in process_names:
                genes = " ".join(df['Target'][df['Source_name'] == process].unique().tolist())
                of.write(process+"\t"+genes+"\n")
                print("WRITTEN PROCESS "+process+"!...")
        print("CREATED GO_DICT FROM HANSEN FORMAT!...")
    
    def jaccardAnalysis(self,t):
        df = t.read('desktop','Mapping.txt','\t',0)
        mbcP = df.MBCO.tolist()
        goP = df.GO.tolist()
        ref_df = t.read('desktop','GO_dict.txt','\t',0)
        data = t.read('downloads','Validation_results_summary_combined.jns','\t',0)
        data.drop(['Description','ReadWrite_identified_terms','ReadWrite_mbc_versions'], axis = 1, inplace = True)
    
        true_positive = ['True_positive','Accepted_positive']
        false_positive = ['False_positive','Denied_positive','Misinterpreted_term','Sibling']
        
        with open('/Users/advaitbalaji/Desktop/JaccardMap.txt','w') as of:
            of.write('MBCGO\tJI\n')
            for process in mbcP:
                genes = data['Symbol'][data['Subcellular_process'] == process].tolist()
                validations = data['Validation'][data['Subcellular_process'] == process].tolist()
                mbcstatus = []
                for v in validations:
                    if v in true_positive:
                        mbcstatus.append('True positive')
                    else:
                        mbcstatus.append('False postitive')
                genes = [gene for gene in genes if mbcstatus[genes.index(gene)] == 'True positive']
                check_process = goP[mbcP.index(process)].split(';')
                for p in check_process:
                    gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                    gene_list[-1] = gene_list[-1].split('\n')[0]
                    of.write(process+'-'+p+'\t'+str(len(list(set(genes) & set(gene_list)))/ len(list(set(genes) | set(gene_list))))+'\n')
      
        data_f = t.read('desktop','JaccardMap.txt','\t',0)
        data_f = data_f[["MBCGO","JI"]][data_f["JI"] != 0]
        ax = sns.boxplot(x="MBCGO", y="JI", data=data_f, whis=np.inf, orient = 'v')
        ax = sns.swarmplot(x="MBCGO", y="JI", data=data_f, color=".2", orient = 'v')
        plt.show()
        
        
                
    def createSummary(self):
        df = t.read('desktop','Mapping.txt','\t',0)
        mbcP = df.MBCO.tolist()
        goP = df.GO.tolist()
        ref_df = t.read('desktop','GO_dict.txt','\t',0)
        data = t.read('downloads','Validation_results_summary_combined.jns','\t',0)
        data.drop(['Description','ReadWrite_identified_terms','ReadWrite_mbc_versions'], axis = 1, inplace = True)
    
        true_positive = ['True_positive','Accepted_positive']
        false_positive = ['False_positive','Denied_positive','Misinterpreted_term','Sibling']
    
        #a,b = t.optionForRIU('both')
    
        with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','w') as metaf:
            metaf.write('MBCProcess\tMissingGenes\n')
    

    
        with open('/Users/advaitbalaji/Desktop/Summary.txt',"w") as sf:
            sf.write("Gene\tMBCO Process\tMBCStatus\tGOStatus\tGO Process description\tIntersection\tUnion\n")
        
            for process in mbcP:
                genes = data['Symbol'][data['Subcellular_process'] == process].tolist()
                validations = data['Validation'][data['Subcellular_process'] == process].tolist()
                #print(goP[mbcP.index(process)])
                check_process = goP[mbcP.index(process)].split(';')
                        
                description = ['' for i in range(len(check_process))]
                GO_megagenelist = []
                #print("This if the first description ",description)
                #print(genes)
                for gene in genes:
                    #print(gene,type(gene))
                    for p in check_process:
                        gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                        gene_list[-1] = gene_list[-1].split('\n')[0]
                        #print(gene in gene_list)
                        GO_megagenelist.extend(gene_list)
                        #if not gene in gene_list:
                            #description[check_process.index(p)] == 'False Positive'
                        if gene in gene_list:
                            description[check_process.index(p)] = 'True Positive'
                        else:
                            description[check_process.index(p)] = 'False Positive'
                    #print(description)
                    #print('')
                    if 'False Positive' in description:
                        if not 'True Positive' in description:
                            gostatus = "Not Present in any mapped GO Process"
                            intersection,union = 'False','False'
                        else:
                            gostatus = "Present in "+str(description.count('True Positive'))+"/"+str(len(description))+" mapped GO Process"
                            intersection,union = 'False','True'
                    else:
                        gostatus = "Present in every mapped GO Process"
                        intersection,union = 'True','True'
                    #print(process,gostatus)              
                    process_descriptions = []
                    for i,j in zip(check_process,description):
                        process_descriptions.append(i+"("+{'False Positive':'Not Present','True Positive': 'Present','': 'False Positive'}[j]+")")
                    process_descriptions = ",".join(process_descriptions)
                    mbcstatus = 'False Positive' if validations[genes.index(gene)] in false_positive else 'True positive'    
                    sf.write(gene+'\t'+process+'\t'+mbcstatus+'\t'+gostatus+'\t'+process_descriptions+'\t'+intersection+'\t'+union+'\n')
            
                with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','a+') as metaf:
                   metaf.write(process+'\t'+','.join(list(set([g for g in GO_megagenelist if not g in genes])))+"\n")    


        



if __name__ == '__main__':
    t = Test()
    lines = t.openFileLineWise('GO_Biological_Process_2017.txt')
    d = t.prepareOrderedDict(lines)  
    d = t.orderDictionaryWithKey(d)
    t.pickleDictToBinary('desktop','pickled_GO',d)
    d = t.readPickledDictFromBinary('desktop','pickled_GO')
    df = t.read('downloads','validation_results_summary_combined.jns','\t',header = 0)
    t.writeDictToFileLineWise(d,'GO_dict')
    
    t.similar('activation of MAPKK activity', 'Acetylcholine-mediated control of postsynaptic potential')
    
    t.createMappingFromBackbone('Downloads/','Mbc_backbone_go_processes_m.txt')
    t.prepareGODictFromJensHansenFormat('Downloads/','gene_association_upgraded_Homo_sapiens_2017June17.jns')
    
    t.jaccardAnalysis(t)
   
