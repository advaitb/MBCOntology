import pandas as pd
from collections import OrderedDict
import pickle
from difflib import SequenceMatcher

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
    
    df = t.read('desktop','Mapping.txt','\t',0)
    mbcP = df.MBCO.tolist()
    goP = df.GO.tolist()
    ref_df = t.read('desktop','GO_dict.txt','\t',0)
    data = t.read('downloads','Validation_results_summary_combined.jns','\t',0)
    data.drop(['Description','ReadWrite_identified_terms','ReadWrite_mbc_versions'], axis = 1, inplace = True)
    
    true_postive = ['True_positive','Accepted_positive']
    false_positive = ['False_positive','Denied_positive','Misinterpreted_term','Sibling']
    
    #a,b = t.optionForRIU('both')
    
    with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','w') as metaf:
        metaf.write('MBCProcess\tMissingGenes\n')
    

    
    with open('/Users/advaitbalaji/Desktop/Summary.txt',"w") as sf:
        sf.write("Gene\tMBCO Process\tMBCStatus\tGOStatus\tGO Process description\tIntersection\tUnion\n")
        
        for process in mbcP:
            genes = data['Symbol'][data['Subcellular_process'] == process].tolist()
            validations = data['Validation'][data['Subcellular_process'] == process].tolist()
            check_process = goP[mbcP.index(process)].split(',')
                        
            description = ['True Positive' for i in range(len(check_process))]
            GO_megagenelist = []
            #print("This if the first description ",description)
            for gene in genes:
                for p in check_process:
                    gene_list = "".join(ref_df['Genes'][ref_df.Process == p].tolist()).split(' ')
                    #print(gene_list)
                    GO_megagenelist.extend(gene_list)
                    #if not gene in gene_list:
                        #description[check_process.index(p)] == 'False Positive'
                    if gene in gene_list:
                        pass
                    else:
                        description[check_process.index(p)] = 'False Positive'
                #print(description)
                if 'False Positive' in description:
                    if not 'True Positive' in description:
                        gostatus = "Not Present in any mapped GO Process"
                        intersection,union = 'False','False'
                    else:
                        gostatus = "Present in "+str(description.count('True Positive'))+"/"+str(len(description))+" mapped GO Process"
                        intersection,union = 'True','False'
                else:
                    gostatus = "Present in every mapped GO Process"
                    intersection,union = 'True','True'
                #print(process,gostatus)              
                process_descriptions = []
                for i,j in zip(check_process,description):
                    process_descriptions.append(i+"("+{'False Positive':'Not Present','True Positive': 'Present'}[j]+")")
                
                process_descriptions = ",".join(process_descriptions)
                mbcstatus = 'False Positive' if validations[genes.index(gene)] in false_positive else 'True positive'    
                sf.write(gene+'\t'+process+'\t'+mbcstatus+'\t'+gostatus+'\t'+process_descriptions+'\t'+intersection+'\t'+union+'\n')
            
            with open('/Users/advaitbalaji/Desktop/missing_genes_per_MBCprocess.txt','a+') as metaf:
               metaf.write(process+'\t'+','.join(list(set([g for g in GO_megagenelist if not g in genes])))+"\n")
                    
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