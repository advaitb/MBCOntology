import xmltodict
from nltk.corpus import stopwords
import nltk
import untangle
#from string import maketrans
#from nltk.stem.snowball import SnowballStemmer
import glob
from sklearn.feature_extraction.text import TfidfVectorizer
import numpy as np
import os
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns
from gensim.corpora import Dictionary
from gensim.models import CoherenceModel, LdaModel, LsiModel, HdpModel
from gensim.models.wrappers import LdaMallet

class NLP:
    
    lemma = nltk.wordnet.WordNetLemmatizer()
    extras = [".",",","=",")","(",";",":","%",">","<","]","[","{","}"]
    path = '/Users/advaitbalaji/Downloads/XML_AbstractNLP/'
    write_path = '/Users/advaitbalaji/Downloads/XML_ProcessedAbstractsNLP/'

    def createDendogramForDengue(self):
        z = np.array([[0.0,0.53,0.22,0.19],[0.53,0.0,0.55,0.63],[0.22,0.55,0.0,0.41],[0.19,0.63,0.41,0.0]])
        df = pd.DataFrame(z,columns =['DENV1','DENV2','DENV3','DENV4'])
        df.set_index([['DENV1','DENV2','DENV3','DENV4']], inplace = True)
        #l_z = linkage(z,"single")
        #data = dendrogram(l_z,orientation = 'top',labels = ['DENV1','DENV2','DENV3','DENV4'], show_leaf_counts = True)
        #print(data)
        #plt.ylabel("HIERARCHIAL CLUSTERING DISTANCE")
        #plt.xlabel("DENGUE VIRUS SEROTYPES")
        #plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6])
        sns.clustermap(df,method = "single",figsize = (5,5), annot = True )
        plt.show()
    
    def is_float(input):
        try:
            num = float(input)
        except ValueError:
            return False
        return True
   
    def tokenize(text):
        return nltk.word_tokenize(text)
    
    def display_scores(vectorizer, tfidf_result):
        # http://stackoverflow.com/questions/16078015/
        #scores = zip(vectorizer.get_feature_names(),
                     #np.asarray(tfidf_result.sum(axis=0)).ravel())
        scores = zip(vectorizer.get_feature_names(),tfidf_result)
        sorted_scores = sorted(scores, key=lambda x: x[1], reverse=True)
        #for item in sorted_scores:
            #return("{0:50} Score: {1}".format(item[0], item[1]))
        sorted_scores = [(a,b) for (a,b) in sorted_scores if b!=0.0 ]
        return('\n'.join([str(a) for a in sorted_scores]))


    def createCorpus(self):
        final_text = ''
        for i in range(243):
            filtered = []
            with open(self.path+'Community_no_'+str(i)+'.txt','r') as f:
                lines = f.readlines()
                for line in lines:
                    if "<AbstractText>" in line:
                        temp_line = line.split('<AbstractText>')[1]
                        temp_line = temp_line.split('</AbstractText>')[0]
                        #abstract = temp_line.translate(trantab)
                        tokens = nltk.word_tokenize(temp_line)
                        filtered_words = [self.lemma.lemmatize(word.lower()) for word in tokens if word not in stopwords.words('english') and word not in self.extras and not is_float(word)]
                        #tokens.remove('.')
                        #tokens.remove(',')
                        filtered = filtered + filtered_words
                final_text = ' '.join(list(set(filtered)))
            with open(self.write_path+'Community_no_'+str(i)+'.txt','w') as w:
                w.write(final_text)
                print("WRITTEN FILE: "+str(i)+"!...")

    def createNLTKHeuristics(self):
        files = sorted(glob.glob(write_path+"*.txt"))
        files = sorted(files, key = lambda x : int(x[69:-4]))

        vector = TfidfVectorizer(input = 'filename', tokenizer = tokenize, analyzer = 'word', stop_words = 'english', ngram_range = (5,5))
        mat = vector.fit_transform(files)
        #print(mat.shape)
        #display_scores(vector,mat)
        df = pd.DataFrame(mat.toarray(), columns = vector.get_feature_names())
        for i in range(243):
            with open('/Users/advaitbalaji/Downloads/Top_hits_NLP(mod)/Community_'+str(i)+'.txt','w') as f:
                f.write(display_scores(vector,df.iloc[i].tolist()))
                print("WRITTEN COMMUNITY : "+str(i)+"!...")
    
                
    def readMultipleFileLineWise(self, file_list):
        texts = []
        for files in file_list:
            with open(files,'r') as f:
                lines = f.readlines()
                texts.append(lines[0].split())
        return texts,file_list
    

    def gensimTopicModelingAnalysis(self,n):
        files  = glob.glob("/Users/advaitbalaji/Downloads/IslandAnalysis/Atleast2/*.txt")
        files = sorted(files, key = lambda x: int(x.split('/Users/advaitbalaji/Downloads/IslandAnalysis/Atleast2/Cluster')[1].split('_')[0]))
        with open("/Users/advaitbalaji/Desktop/ListofSortedClusters.txt","w") as of:
            for f in files:
                of.writelines(f+"\n")
        texts,clusters = n.readMultipleFileLineWise(files)
        dictionary = Dictionary(texts)
        corpus = [dictionary.doc2bow(text) for text in texts]
        hdpmodel = HdpModel(corpus=corpus, id2word=dictionary)
        print(hdpmodel.show_topics())
            #print()

        





        
if __name__ == '__main__':
    n = NLP()
    n.gensimTopicModelingAnalysis(n)
    





'''
p = doc['PubmedArticleSet']['Article']

info = list(p[0].items())
a,b = info[len(info)-1]
print(b)
print('\n')
print('\n')
words = nltk.word_tokenize(b)
for word in words:
    if word in stopwords.words('english'):
        words.remove(word)
#words.remove('.').remove(',')
freq = nltk.FreqDist(words)
for key,val in freq.items():
    print(str(key) + ':' + str(val))

#print(words)

''' 
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns
#plt.style.use('ggplot')

z = np.array([[0.0,0.53,0.22,0.19],[0.53,0.0,0.55,0.63],[0.22,0.55,0.0,0.41],[0.19,0.63,0.41,0.0]])
df = pd.DataFrame(z,columns =['DENV1','DENV2','DENV3','DENV4'])
df.set_index([['DENV1','DENV2','DENV3','DENV4']], inplace = True)
#l_z = linkage(z,"single")

#data = dendrogram(l_z,orientation = 'top',labels = ['DENV1','DENV2','DENV3','DENV4'], show_leaf_counts = True)
#print(data)
#plt.ylabel("HIERARCHIAL CLUSTERING DISTANCE")
#plt.xlabel("DENGUE VIRUS SEROTYPES")
#plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6])

sns.clustermap(df,method = "single",figsize = (5,5), annot = True )
plt.show()
'''
