from Bio import Entrez
import time
import datetime as dt
import os
import pandas as pd
def search(query):
    Entrez.email = 'advait.balaji@mssm.edu'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='999999',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list): 
    ids = ','.join(id_list)
    Entrez.email = 'advait.balaji@mssm.edu'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

#il_list = [str(j) for j in range(1,14)]
#for i in il_list:
with open('/Users/advaitbalaji/Desktop/untitled.txt','r') as f:
    
    lines = f.readlines()
    
    for line in lines:
        
        print(line)
        
        if line == '\n':
            continue
        
        fname = line.replace('-','').split('\t')[0]
        searchq = line.replace('-','').split('\t')[1]
    
        with open('/Users/advaitbalaji/Downloads/NewProcesses2018/XML_'+fname+'.txt','w') as of:
            
                results = search(searchq)
                id_list = results['IdList']
                print(searchq)
                
                if len(id_list) == 0:
                    continue
                if len(id_list) > 5000:
                    id_list = id_list[:5000]

                #id_list = [str(a) for a in df[df['modularity_class'] == i]['Id'].tolist()]
                papers = fetch_details(id_list)
                of.write('<PubmedSubcellularProcessDownload>\n')
                of.write('\t<General_Information>\n')
                of.write('\t\t<subcellular_process_name>'+fname+'</subcellular_process_name>\n')
                of.write('\t\t<Pubmed_querry>'+'Abstract_Check_'+searchq+'</Pubmed_querry>\n')
                of.write('\t\t<Pubmed_article_count>'+str(len(id_list))+'</Pubmed_article_count>\n')
                of.write('\t\t<Download_date>' + dt.datetime.today().strftime("%m/%d/%Y") + '</Download_date>\n')
                of.write('\t</General_Information>\n')
                of.write('\t<PubmedArticleSet>\n')
            	            #print(papers)s
                for paper in papers['PubmedArticle']:
                    if 'PMID' not in paper['MedlineCitation']:
                        pmid = 'NOT AVAILABLE'
                    else:
                        pmid = paper['MedlineCitation']['PMID']
                    if 'ArticleTitle' not in paper['MedlineCitation']['Article']:
                        article_title = 'NOT AVAILABLE'
                    else:
                        article_title = paper['MedlineCitation']['Article']['ArticleTitle']
                    journal= paper['MedlineCitation']['Article']['Journal']['Title']
                    if 'PubDate' not in paper['MedlineCitation']['Article']['Journal']['JournalIssue']:
                        pubdate = 'NOT AVAILABLE'
                    else:
                        if 'Year' not in paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                            pubdate = 'NOT AVAILABLE'
                        elif 'Month' not in paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                            pubdate = 'NOT AVAILABLE'
                        else:
                            pubdate = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year'] + '-' + paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
                    if 'Abstract' in paper['MedlineCitation']['Article']:
                        abstract = ''.join(paper['MedlineCitation']['Article']['Abstract']['AbstractText'])
                    else:
                        abstract = 'NOT AVAILABLE'
                    of.write('\t\t<Article>\n')
                    of.write('\t\t\t<PMID>'+pmid+ '</PMID>\n')
                    of.write('\t\t\t<title>'+article_title+ '</title>\n')
                    of.write('\t\t\t<JournalAbbreviation>'+journal+'</JournalAbbreviation>\n')
                    of.write('\t\t\t<Publication_date>'+pubdate+ '</Publication_date>\n')
                    of.write('\t\t\t<AbstractText>'+abstract+'</AbstractText>\n')
                    of.write('\t\t</Article>\n')
                of.write('\t</PubmedArticleSet>\n')
                of.write('</PubmedSubcellularProcessDownload>')
                print("Written Process: "+ fname + "....")
    
	#print('Toll_like_receptor '+i+' signalling XML file done...')


    #for i, paper in enumerate(papers['PubmedArticle']):
    #   print("%d) %s" % (i+1, paper['MedlineCitation']['Article']['ArticleTitle']))
    # Pretty print the first paper in full
    #import json
    #print(json.dumps(papers[0], indent=2, separators=(',', ':')))
