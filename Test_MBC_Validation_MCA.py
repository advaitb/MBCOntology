import pandas as pd
from collections import OrderedDict
import pickle
from difflib import SequenceMatcher
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances,silhouette_score,calinski_harabaz_score

class Test:
    
    path = '/Users/advaitbalaji/'
    downloads = 'Downloads/'
    desktop = 'Desktop/'
    celltypes_group = ['Vascular endothelial cell','T cell','Stromal cell','Stomach cell','Stem and progenitor cell','Spongiotrophoblast','Spermatogonia','Spermatocyte','Spermatids','Smooth muscle cell','Skeletal muscle cell','Sertoli cell','Secretory alveoli cell','Schwann cell','S cell','Radial glia','Pyramidal neuron cell','Pre-Sertoli cell','Pit cell','PE lineage cell','Osteoblast','Neutrophil','Neuron','NK cell','Myoblast','Muscle  cell','Monocyte','Microglia','Mesenchymal stem cell','Mesenchymal Alveolar Niche Cell','Megakaryocyte','Mast cell','Marcrophage','Luminal cell','Leydig cell','Keratinocyte','Hypothalamic ependymal cell','Hepatocyte','Granulosa cell','Granulocyte','Glomerular','Glial cell','Glandular epithelium','Ganglion cell','Fibroblasts','Erythrocyte','Erythroblast','Epithelium of small intestinal villi','Epithelial cell','Eosinophil','Enterocyte progenitor','Endothelial cells','Endocrine cell','ES_Rps28_high','Dividng cell','Dendrtic cell','Cumulus cell','Conventional dendritic cell','Cartilage cell','Brown adipose tissue','Bipolar cell','Basophil','B cell','Axin2+ Myofibrogenic Progenitor cell','Astrocyte','Amacrine cell','Alveolar macrophage','Adipocyte','Acinar cell']
    
    def checkMBCSCP(self,inp):
        #TAB SEPARATED FILE AS INPUT, TARGET IS THE PROCESS WHERE FALSE_POSTITVE WILL BE ASSIGNED!..
        df = pd.read_csv('/Users/advaitbalaji/Desktop/Validation_results_summary_combined.jns', sep = '\t', header = 0)
        with open(self.path+inp,'r') as inpf:
            print("Correcting MBC Validation summary!...")
            lines = inpf.readlines()
            for line in lines:
                source, target = line.split('\t')
                s_genes = df[(df['Subcellular_process'] == source) & (df['Validation'] == 'True_positive')]['Symbol'].tolist()
                t_genes = df[df['Subcellular_process'] == source]['Symbol'].tolist()
                check_genes = [g for g in s_genes if g in t_genes]
                for gene in check_genes:
                    df.loc[(df['Symbol'] == gene) & (df['Subcellular_process'] == target),'Validation'] = 'False_positive'
                df.to_csv('/Users/advaitbalaji/Desktop/Validation_results_summary_combined_CORRECTED.jns', sep = '\t', header = True, index = False)
    
    def makeListOfCellTypes(self):
        #LEVEL==1 FOR GO AND LEVEL==3 FOR MBC
        celltypes = {}
        processes = []
        all_columns = ['Processes']
        files = glob.glob('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/*.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/Arc-Me(Drop-seq)_SingleCell.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/E18-Brain(10X)_SingleCell.csv')
        for f in files:
            key = f.split('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/')[1].split('_SingleCell.csv')[0]
            #print(key)    
            df = pd.read_csv(f, sep = '\t', header = 0)
            #celltypes.extend([cell+':'+key for cell in df.SampleName.unique().tolist()])
            celltypes[key] = [cell for cell in df.SampleName.unique().tolist()]
        for k in celltypes.keys():
            
            for cells in celltypes[k]:
                #print(type(cells))
                all_columns.extend([cells+':'+k])
            data = pd.read_csv('/Users/advaitbalaji/Desktop/MBCO_enrichment_current/Results/'+k+'_SingleCell.csv_Molecular_biology_cell/Standard_enrichment_results_filtered.txt', sep = '\t', header = 0)
            #USE BOTTOM 2 ONLY FOR GO ANALYSIS NOT FOR MBCO!...
            data = data[['ProcessLevel','Sample_name','Scp_name','Minus_log10_pvalue']]
            data = data.groupby(['Sample_name']).apply(lambda x: x.sort_values(['Minus_log10_pvalue'], ascending = False).head(10)).reset_index(drop=True)
            processes.extend(data[data['ProcessLevel'] == 1]['Scp_name'].unique().tolist())
        processes = list(set(processes))
        #print(len(all_columns))
        #print(processes)
        with open('/Users/advaitbalaji/Desktop/go_mcaanalysis.txt','a') as of:
            for process in processes:
                of.write(process+''.join(['\t' for i in range(len(all_columns)-1)])+'\n')
            
        #with open('/Users/advaitbalaji/Desktop/mcaanalysis.txt','w') as of:
        #    of.write('\t'.join(all_columns)+'\n')
        #SORT DATAFRAME ACCORDING TO CELL TYPE AND KEEP MAKING DICT OF THE TOP 10 PROCESSES, IF LESS THAN 10 THEN TAKE WHATEVER AVAILABLE!...
    
    #GENE LIST IS NONSENSE!
    def makeGenesList(self):
        celltypes = {}
        genes = []
        all_columns = ['Genes']
        files = glob.glob('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/*.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/Arc-Me(Drop-seq)_SingleCell.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/E18-Brain(10X)_SingleCell.csv')
        for f in files:
            key = f.split('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/')[1].split('_SingleCell.csv')[0]
            #print(key)    
            df = pd.read_csv(f, sep = '\t', header = 0)
            #celltypes.extend([cell+':'+key for cell in df.SampleName.unique().tolist()])
            celltypes[key] = [cell for cell in df.SampleName.unique().tolist()]
        for k in celltypes.keys():
            for cells in celltypes[k]:
                #print(type(cells))
                all_columns.extend([cells+':'+k])
            data = pd.read_csv('/Users/advaitbalaji/Desktop/MBCO_enrichment_current/Results/'+k+'_SingleCell.csv_Molecular_biology_cell/Standard_enrichment_results_filtered.txt', sep = '\t', header = 0)
            #data = data[['ProcessLevel','Sample_name','Scp_name','Minus_log10_pvalue']]
            #data = data.groupby(['Sample_name']).apply(lambda x: x.sort_values(['Minus_log10_pvalue'], ascending = False).head(10)).reset_index(drop=True)
            genes.extend([b.split(' ')[0] for a in data[data['ProcessLevel'] == 1]['ReadWrite_overlap_symbols'].tolist() for b in a.split(',')])
        genes = list(set(genes))
        #print(len(all_columns))
        #print(processes)
        with open('/Users/advaitbalaji/Desktop/go_mcaanalysis-genes.txt','a') as of:
            for gene in genes:
                of.write(gene+''.join(['\t' for i in range(len(all_columns)-1)])+'\n')
                
   #GENE CLUSTER IS NONSENSE!
    def setGeneClusterAnalysisFile(self,f):
        df = pd.read_csv(f,sep = '\t',header = 0)
        #processes = df.Processes.tolist()
        genes = df.Genes.tolist()
        columns = list(df.columns)[1:]
        for column in columns:
            cell,tissue = column.split(':')
            #print(type(cell),type(tissue))
            data = pd.read_csv('/Users/advaitbalaji/Desktop/MBCO_enrichment_current/Results/'+tissue+'_SingleCell.csv_Molecular_biology_cell/Standard_enrichment_results_filtered.txt', sep = '\t', header = 0)
            #BELOW ONLY FOR GO!...
            #data = data[['ProcessLevel','Sample_name','Scp_name','Minus_log10_pvalue']]
            #data = data.groupby(['Sample_name']).apply(lambda x: x.sort_values(['Minus_log10_pvalue'], ascending = False).head(10)).reset_index(drop=True)
            g = [b.split(' ')[0] for a in data[(data['ProcessLevel'] == 1) & (data['Sample_name'] == cell)]['ReadWrite_overlap_symbols'].tolist() for b in a.split(',')]
            pvalue = list(-np.log10([float(b.split(' ')[1].strip(')').strip('(')) for a in data[(data['ProcessLevel'] == 1) & (data['Sample_name'] == cell)]['ReadWrite_overlap_symbols'].tolist() for b in a.split(',')]))
            #print(pvalue)
            #break
            for gene in genes:
                if gene in g:
                     df.loc[df['Genes'] == gene,column] = pvalue[g.index(gene)]
                else:
                    df.loc[df['Genes'] == gene,column] = 0.0
            print('WRITTEN CELL TYPE: '+str(cell)+' IN TISSUE: '+str(tissue)+' '+str(columns.index(column)+1)+'/'+str(len(columns))+'!...')
            
        df.to_csv(f,sep='\t',header=True,index=False)
     
    
    def setClusterAnalysisFile(self,f):
        df = pd.read_csv(f,sep = '\t',header = 0)
        #processes = df.Processes.tolist()
        processes = df.Processes.tolist()
        columns = list(df.columns)[1:]
        for column in columns:
            cell,tissue = column.split(':')
            #print(type(cell),type(tissue))
            data = pd.read_csv('/Users/advaitbalaji/Desktop/MBCO_enrichment_current/Results/'+tissue+'_SingleCell.csv_Molecular_biology_cell/Standard_enrichment_results_filtered.txt', sep = '\t', header = 0)
            #BELOW ONLY FOR GO!...
            #data = data[['ProcessLevel','Sample_name','Scp_name','Minus_log10_pvalue']]
            #data = data.groupby(['Sample_name']).apply(lambda x: x.sort_values(['Minus_log10_pvalue'], ascending = False).head(10)).reset_index(drop=True)
            scps = data[(data['ProcessLevel'] == 1) & (data['Sample_name'] == cell)]['Scp_name'].tolist()
            pvalue = data[(data['ProcessLevel'] == 1) & (data['Sample_name'] == cell)]['Minus_log10_pvalue'].tolist()
            
            for process in processes:
                if process in scps:
                     df.loc[df['Processes'] == process,column] = pvalue[scps.index(process)]
                else:
                    df.loc[df['Processes'] == process,column] = 0.0
            print('WRITTEN CELL TYPE: '+str(cell)+' IN TISSUE: '+str(tissue)+' '+str(columns.index(column)+1)+'/'+str(len(columns))+'!...')
            
        df.to_csv(f,sep='\t',header=True,index=False)
        #print(list(df.columns))
    
    def makeClusterMap(self,f):
        df = pd.read_csv(f,sep = '\t', header = 0, index_col = 'Processes')
        #print(len(df.columns))
        #print(len(df.values[0]))
        ac = AgglomerativeClustering(n_clusters = 100)
        clusters = ac.fit_predict(df.values.T)
        print(clusters)
        cluster_dict = {i:[] for i in range(100)}
        for key,column in enumerate(df.columns):
            cluster_dict[clusters[key]].append(column)
            #print(key)
        print(cluster_dict)
        #plt.scatter(clusters,ac.n_clusters)
        #plt.show()
        #print(df.values),ac.n_clusters
        #pca = PCA()
        #pca.fit_transform(df.values)
        #sns.set(font_scale=0.5)
        #g = sns.clustermap(df, cmap = 'Blues', figsize = (13,13))
        #plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)  # For y axis
        #plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90) # For x axis
        #plt.xticks(df.columns[1:],rotation=90)
        #plt.yticks(df.Processes.tolist(),rotation = 90)
        #plt.tick_params(labelsize=2)
        #g.ax_heatmap.tight_layout()
        #print(g.dendrogram_col.linkage)
        #plt.show()
        
    def findPairWiseDistance(self):
        df = pd.read_csv('/Users/advaitbalaji/Desktop/genesmca.txt',sep = '\t', header = 0, index_col = 'Genes')
        #with open(self.path+self.desktop+'CellTypes.txt','w') as of:
        #    for column in sorted(df.columns[1:]):
        #        of.write(column+'\n')
        #print(df.isnull().sum().sum())
        values = df.values.T
        intracelldist = []
        intercelldist = []
        cell_index = []
        for index,cell in enumerate(df.columns):
            if 'Granulocyte' in cell:
                #print(cell)
                cell_index.append(index)
        #print(len(cell_index))
        other_index = [index for index in range(803) if index not in cell_index]
        #print(len(other_index))
        #print(cell_index)
        for i in range(len(cell_index)):
            for j in range(i+1,len(cell_index)):
                value1 = np.reshape(values[i],(1,-1))
                value2 = np.reshape(values[j],(1,-1))
                intracelldist.append(pairwise_distances(value1,value2, metric = "euclidean").mean())
                #intracelldist.append(np.nanmean(pairwise_distances(value1,value2, metric = "euclidean")))
        #print(intracelldist)
        #print(np.mean(intracelldist))
        for c in cell_index:
            for o in other_index:
                value1 = np.reshape(values[c],(1,-1))
                value2 = np.reshape(values[o],(1,-1))
                intercelldist.append(pairwise_distances(value1,value2, metric = "euclidean").mean())
                #intercelldist.append(np.nanmean(pairwise_distances(value1,value2, metric = "euclidean")))
        #print(intercelldist)
        print((np.mean(intercelldist) - np.mean(intracelldist))/max(np.mean(intercelldist),np.mean(intracelldist)))
              
        '''
        data = pd.read_csv('/Users/advaitbalaji/Desktop/go_mcaanalysis.txt')
        values = df.values.T
        #WE CAN ITERATE OVER THE NEEDED CELLS BY REFERRING TO THE INDICES IN DF.COLUMNS, CALCULATE PW DISTANCE, ALSO INTER-CLUSTER AND INTRA CLUSTER DISTANCE!..
        value1 = np.reshape(values[1],(1,-1))
        value2 = np.reshape(values[200],(1,-1))
        distance = pairwise_distances(value1,value2, metric = "euclidean").mean() 
        print(distance)
        '''
        
    def geneBasedClustering(self):
        celltypes = {}
        genes = []
        all_columns = ['Genes']
        files = glob.glob('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/*.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/Arc-Me(Drop-seq)_SingleCell.csv')
        files.remove('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/E18-Brain(10X)_SingleCell.csv')
        for f in files:
            key = f.split('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/')[1].split('_SingleCell.csv')[0]
            #print(key)    
            df = pd.read_csv(f, sep = '\t', header = 0)
            #celltypes.extend([cell+':'+key for cell in df.SampleName.unique().tolist()])
            celltypes[key] = df.SampleName.unique().tolist()
            genes.extend(df.NCBI_official_symbol.unique().tolist())
        genes = list(set(genes))    
        for k in celltypes.keys():
            for cells in celltypes[k]:
                #print(type(cells))
                all_columns.extend([cells+':'+k])
        with open('/Users/advaitbalaji/Desktop/genesmca.txt','w') as of:
            for gene in genes:
                of.write(gene+''.join(['\t' for i in range(len(all_columns)-1)])+'\n')
    
    def setGenesMCA(self,f):
        df = pd.read_csv(f,sep = '\t',header = 0)
        genes = df.Genes.tolist()
        columns = list(df.columns)[1:]
        for column in columns:
            cell,tissue = column.split(':')
            data = pd.read_csv('/Users/advaitbalaji/Desktop/MCA_EnrichmentFormat/'+tissue+'_SingleCell.csv', sep = '\t', header = 0)
            g = data[data['SampleName'] == cell]['NCBI_official_symbol'].tolist()
            pvalue = list(-np.log10([float(a) if float(a) != 0.0 else 1.0 for a in data[data['SampleName'] == cell]['Value'].tolist()]))
            
            for gene in g:
                df.loc[df['Genes'] == gene,column] = pvalue[g.index(gene)]
            #else:
            #        df.loc[df['Genes'] == gene,column] = 0.0
                    
            print('WRITTEN CELL TYPE: '+str(cell)+' IN TISSUE: '+str(tissue)+' '+str(columns.index(column)+1)+'/'+str(len(columns))+'!...')
        df.replace([np.inf, -np.inf], np.nan, inplace = True)
        df.fillna(value = 0.0, inplace = True)
        df.to_csv(f,sep='\t',header=True,index=False) 
        
    
    
    def findSilhouetteScore(self,f):
        df = pd.read_csv(f,sep = '\t', header = 0, index_col = 'Processes')
        ac = AgglomerativeClustering(n_clusters = 100)
        clusters = ac.fit_predict(df.values.T)
        score = silhouette_score(df.values.T,clusters,metric='euclidean')
        print(score)
    
    def findClinskiharabazScore(self,f):
        df = pd.read_csv(f,sep = '\t', header = 0, index_col = 'Processes')
        ac = AgglomerativeClustering(n_clusters = 100)
        clusters = ac.fit_predict(df.values.T)
        score = calinski_harabaz_score(df.values.T,clusters)
        print(score)
        #calinski_harabaz_score 
    
    
    def makeGO_dictSuppReady(self,f):
        df = pd.read_csv(f, sep = '\t', header = 0)
        with open('/Users/advaitbalaji/Desktop/Supplementary Table S32 - gene-SCP associations.txt','w+') as of:
            of.write("ProcessLevel\tParent_processName\tProcessID\tProcessName\tSymbol\n")
            for process in df.Process.tolist():
                g = df['Genes'][df['Process'] == process].tolist()
                genes = g[0].split(' ')
                for gene in genes:
                    of.write('1\tGO_Parent\tGO_PID\t'+process+'\t'+gene+'\n')
    
    
        
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
                    of.write(process+'::'+p+'\t'+str(len(list(set(genes) & set(gene_list)))/ len(list(set(genes) | set(gene_list))))+'\n')
        
        
        #t.prepareLevel21correlation()
        #t.prepareLevel32correlation()
        
    def createLevel1ProcessMappedJI(self,t):
        df_32 = t.read('desktop','Level32.csv','\t',0)
        df_21 = t.read('desktop','Level21.csv','\t',0)
            
        with open('/Users/advaitbalaji/Desktop/JaccardMap.txt','r') as inpf,open('/Users/advaitbalaji/Desktop/JaccardMapParent.txt','w') as outpf:
            outpf.write('Level1ProcessMBCO\tJI\n')
            lines = inpf.readlines()
            for line in lines[1:]:
                a,b = line.split('\t')
                process = a.split('::')[0]
                #print(process)
                level1p = "".join(df_21['Parent_process'][df_21['Process'] == "".join(df_32['Parent_process'][df_32['Process'] == process].tolist()) ].tolist())
                if level1p == "":
                    level1p = "".join(df_21['Parent_process'][df_21['Process'] == process].tolist())
            
                #print(str(lines.index(line))+". "+level1p)
                #print("")
                if not level1p == "":
                    outpf.write(level1p+'\t'+b)
                else:
                    continue
                
        #data_f = t.read('desktop','JaccardMap.txt','\t',0)
        #data_f = data_f[["MBCGO","JI"]][data_f["JI"] != 0]
        #ax = sns.boxplot(x="MBCGO", y="JI", data=data_f, whis=np.inf, orient = 'v')
        #ax = sns.swarmplot(x="MBCGO", y="JI", data=data_f, color=".2", orient = 'v')
        #plt.show()
    
    def createBoxSwarmPlot(self,inputf):
        data = pd.read_csv(self.path+self.desktop+inputf, sep = '\t', header = 0)
        plt.figure(figsize=(7,7))
        ax = sns.boxplot(x="Level1ProcessMBCO", y="JI", data=data, whis=np.inf)
        ax = sns.swarmplot(x="Level1ProcessMBCO", y="JI", data=data, color=".2")
        #plt.figure(figsize=(20,10))
        plt.xlabel("Level 1 Processes in MBC Ontology")
        plt.ylabel("Jaccard Index")
        plt.xticks(rotation=90)
        plt.tick_params(axis='x', which='major', labelsize=7)
        plt.tight_layout()
        plt.show()
        
    def prepareLevel21correlation(self):
        data = pd.read_csv('/Users/advaitbalaji/Desktop/Mbc_backbone_regular.csv', sep = '\t', header = 0)
        data = data[['Process','Parent_process']][data['Level'] == 2]
        data = data.reset_index(drop=True)
        data.to_csv('/Users/advaitbalaji/Desktop/Level21.csv', sep = '\t', header = True, index  = False)
    
    def prepareLevel32correlation(self):
        data = pd.read_csv('/Users/advaitbalaji/Desktop/Mbc_backbone_regular.csv', sep = '\t', header = 0)
        data = data[['Process','Parent_process']][data['Level'] == 3]
        data = data.reset_index(drop=True)
        data.to_csv('/Users/advaitbalaji/Desktop/Level32.csv', sep = '\t', header = True, index  = False)
            
                
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
    #lines = t.openFileLineWise('GO_Biological_Process_2017.txt')
    #d = t.prepareOrderedDict(lines)  
    #d = t.orderDictionaryWithKey(d)
    #t.pickleDictToBinary('desktop','pickled_GO',d)
    #d = t.readPickledDictFromBinary('desktop','pickled_GO')
    #df = t.read('downloads','validation_results_summary_combined.jns','\t',header = 0)
    #t.writeDictToFileLineWise(d,'GO_dict')
    
    #t.similar('activation of MAPKK activity', 'Acetylcholine-mediated control of postsynaptic potential')
    
    #t.createMappingFromBackbone('Downloads/','Mbc_backbone_go_processes_m.txt')
    #t.prepareGODictFromJensHansenFormat('Downloads/','gene_association_upgraded_Homo_sapiens_2017June17.jns')
    
    #t.jaccardAnalysis(t)
    #t.createLevel1ProcessMappedJI(t)
    #t.createBoxSwarmPlot('JaccardmapParent.txt')
    #t.checkMBCSCP('Desktop/Input.txt')
    #t.makeListOfCellTypes()
    #t.setClusterAnalysisFile('/Users/advaitbalaji/Desktop/go_mcaanalysis.txt') 
    #t.makeClusterMap('/Users/advaitbalaji/Desktop/go_mcaanalysis.txt')
    #t.findSilhouetteScore('/Users/advaitbalaji/Desktop/mcaanalysis.txt')
    #t.findClinskiharabazScore('/Users/advaitbalaji/Desktop/go_mcaanalysis.txt')
    t.findPairWiseDistance()
    #t.geneBasedClustering()
    #t.setGenesMCA('/Users/advaitbalaji/Desktop/genesmca.txt')
    #t.makeGenesList()
    #t.setGeneClusterAnalysisFile('/Users/advaitbalaji/Desktop/go_mcaanalysis-genes.txt')       
    #t.makeGO_dictSuppReady('/Users/advaitbalaji/Downloads/MBCMapping/GO_dict.txt')

'''
#Frequent Words Function
def FrequentWordsWithMismatches( s, k, d ):
    counts = {}
    for i in range(len(s)-k+1):
        for neighbor in neighbors(s[i:i+k],d):
            print(neighbor)
            if neighbor not in counts:
                counts[neighbor] = 0
            counts[neighbor] += 1
    m = max(counts.values())
    return [kmer for kmer in counts if counts[kmer] == m]

#Finding neighbouring kmers
def neighbors( s, d ):
    #print(s, d)
    if d == 0:
        return [s]
    if len(s) == 1:
        return ['A','C','G','T']
    out = []
    for neighbor in neighbors(s[1:],d):
        if hamming(s[1:],neighbor) < d:
            out.extend(['A'+neighbor,'C'+neighbor,'G'+neighbor,'T'+neighbor])
        else:
            out.append(s[0] + neighbor)
    print(out)    
    return out

#Hamming Distance between kmers
def hamming( s, t ):
    return sum([s[i] != t[i] for i in range(len(s))])

#s = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
s = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
k = 4
d = 1
print(' '.join(FrequentWordsWithMismatches(s,k,d)))

'''



'''
sequence = 'CATTCCAGTACTTCATGATGGCGTGAAGA'

skew = []
prev_val = 0
skew.append(prev_val)
for i in sequence:
    if i == 'C':
        prev_val -= 1
        skew.append(prev_val)
    elif i == 'G':
        prev_val += 1
        skew.append(prev_val)
    else:
        skew.append(prev_val)
#print(skew)
indexes = [i for i,j in enumerate(skew) if j == max(skew)]
print(indexes)
'''

'''
prev_val = 0
min_val = 0
indexes = []
indexes.append(0)
for i in range(len(Genome)):
    if Genome[i] == 'C':
        prev_val -= 1
        if prev_val < min_val:
            indexes = []
            indexes.append(i+1)
            min_val = prev_val
        elif prev_val == min_val:
            indexes.append(i+1)
        else:
            pass            
    elif Genome[i] == 'G':
        prev_val += 1
    else:
        pass
print(indexes)
'''







'''
with open("/Users/advaitbalaji/Desktop/EcGenome.txt","r") as f:
    genome = f.read()
    
k, Length , t = 9,500,3
overall_count = []
for i in range(len(genome) - Length + 1):
    count_dict ={}
    mod_genome = genome[i:i+Length]
    for j in range(len(mod_genome) - k + 1):
        if mod_genome[j:j+k] in count_dict:
            count_dict[mod_genome[j:j+k]] += 1
            if count_dict[mod_genome[j:j+k]] >= t:
                overall_count.append(mod_genome[j:j+k])
        else:
            count_dict[mod_genome[j:j+k]] = 1
#print("YES")
print(" ".join(list(set(overall_count))))
'''            


'''
with open("/Users/advaitbalaji/Desktop/VBGenome.txt","r") as f:
    Genome = f.read()

Pattern = "CTTGATCAT"



pat_index = []
for i in range(len(Genome) - len(Pattern) + 1):
    if Genome[i:i+len(Pattern)] == Pattern:
        pat_index.append(str(i))
print(" ".join(pat_index))

'''

'''
Text = 'CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA'
for k in range(3,4):
    count_dict = {}
    max_count = 0
    for i in range(len(Text) - k + 1):
        if Text[i:i+k] in count_dict:
            count_dict[Text[i:i+k]] += 1
            if max_count < count_dict[Text[i:i+k]]:
                max_count = count_dict[Text[i:i+k]]
        else:
            count_dict[Text[i:i+k]] = 1

    max_kmers = []
    for i, j in count_dict.items():
        if j == max_count:
            max_kmers.append(i)
    print(str(k)+": "+" ".join(max_kmers)+ " "+ str(max_count))
'''