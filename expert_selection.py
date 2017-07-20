""" Fetch author emails from Pubmed
:Author: Jonathan Karr <karr@mssm.edu> & Balazs Szigeti
:Date: 2017-07-12
:Copyright: 2016, Karr Lab
:License: MIT
"""

from Bio import Entrez
from difflib import SequenceMatcher

import csv
import datetime
import itertools
import math
import numpy
import scholarly
import time

# configure entrez
Entrez.email = 'balazs.szigeti@mssm.edu'
max_results = 1000
query = '("Cell Biology/analysis"[Mesh] OR "Cell Biology/drug effects"[Mesh] OR "Cell Biology/genetics"[Mesh] OR "Cell Biology/physiology"[Mesh] OR "Cell Biology/statistics and numerical data"[Mesh]) AND "Computational Biology"[Majr] OR "Computer Simulation"[Majr]'

print('\n','Query: \n',query,'\n')

def main():
    # search PubMed
    n_error        = 0
    search_handle  = Entrez.esearch(db='pubmed', sort='relevance', retmax=str(max_results), term=query)
    search_results = Entrez.read(search_handle)
    pubmed_ids     = search_results['IdList']

    # Fetch publication metadata
    corresponding_authors = []
    start_time = time.time()
    chunk_size = 25

    for i_chunk in range(int(math.ceil(len(pubmed_ids) / chunk_size))):
        print('Retrieving articles {}..{} of {}'.format(i_chunk * chunk_size + 1, min(len(pubmed_ids), (i_chunk + 1) * chunk_size), len(pubmed_ids)))

        fetch_handle        = Entrez.efetch(db='pubmed', retmode='xml', id=','.join(pubmed_ids[i_chunk * chunk_size:min(len(pubmed_ids), (i_chunk + 1) * chunk_size)]))
        entrez_query_return = Entrez.read(fetch_handle)

        for article in entrez_query_return['PubmedArticle']:
            Pubmed_title    = article['MedlineCitation']['Article']['ArticleTitle']
            GS_search_query = scholarly.search_pubs_query(Pubmed_title )

            try:
                GS_article      = next(GS_search_query)
            except:
                print('weird error')
                continue

            GS_title        = GS_article.get_title()

            title_distance = SequenceMatcher(None,Pubmed_title, GS_title).ratio()
            if title_distance > 0.8:
                corresponding_authors.extend(parse_entrez_query_return(article,GS_article,n_error))
            else:
                n_error = n_error+1
                print('Unmatching titles, title distance: ',round(title_distance,3))

    # Group authors by email
    print('done')
    investigators = []
    key = lambda author: author['email'].lower()
    corresponding_authors.sort(key=key)


    for email, articles in itertools.groupby(corresponding_authors, key):
        sorted_articles = sorted(articles,key=lambda article: article['publication_date'] if article['publication_date'] else datetime.date(1, 1, 1),reverse=True)
        lastest_article = sorted_articles[0]

        investigators.append({
            'first_name': lastest_article['first_name'],
            'last_name': lastest_article['last_name'],
            'No of articles': len(sorted_articles),
            'citations': [article['n_citation'] for article in sorted_articles],
            'abstract': lastest_article['abstract'],
            'email': email,
            'pubmed_ids': [article['pubmed_id'] for article in sorted_articles],
            #'affiliation': lastest_article['affiliation'],
            #'all_affiliations': [article['affiliation'] for article in sorted_articles],
            #'lastest_publication_date': lastest_article['publication_date'],
            })

    # sort investigators by publication number
    #investigators.sort(key=lambda investigator: (investigator['No of articles'], investigator['last_name'], investigator['first_name']),reverse=True)
    investigators.sort(key=lambda investigator: (investigator['citations'], investigator['No of articles'], investigator['last_name']),reverse=True)

    # count number of authors
    article_count_frequencies = numpy.bincount([len(i['pubmed_ids']) for i in investigators])
    print('\n')
    print('{:12}  {:17}'.format('No. articles', 'No. investigators'))
    for n_articles, n_investigators in enumerate(article_count_frequencies):
        print('{:12}  {:17}'.format(n_articles, n_investigators))

    print(str(len(pubmed_ids))+' papers were retrived, out of which '+str(n_error)+' could not be processed.')

    # save author information to a tsv file
    with open('investigators.tsv', 'w+') as file:
        csv_writer = csv.DictWriter(file, delimiter='\t', fieldnames=investigators[0].keys())
        csv_writer.writeheader()
        for investigator in investigators:
            csv_writer.writerow(investigator)

        file.close()
        print("--- %s seconds ---" % (time.time() - start_time))

def parse_entrez_query_return(article,GS_article,n_error):

    try:
        abstract_text = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
    except:
        n_error = n_error+1
        return []

    medline_article = article['MedlineCitation']['Article']
    if 'AuthorList' not in medline_article:
        n_error = n_error+1
        return []

    author_list = medline_article['AuthorList']
    if not author_list:
        n_error = n_error+1
        return []

    corresponding_authors = []
    for author in author_list:
        affiliation_info = author['AffiliationInfo']
        if not affiliation_info:
            n_error = n_error+1
            continue

        affiliation = affiliation_info[0]['Affiliation']
        affiliation, _, email = affiliation[0:-1].rpartition(' ')
        if '@' not in email:
            n_error = n_error+1
            continue

        id_list = article['PubmedData']['ArticleIdList']
        pubmed_id = None
        for id in id_list:
            if id.attributes['IdType'] == 'pubmed':
                pubmed_id = str(id)

        try:
            n_citation = GS_article.get_citation()
        except:
            n_citation = 0 # If there are no citation, they are not founs

        if 'ArticleDate' in medline_article and medline_article['ArticleDate']:
            tmp = medline_article['ArticleDate'][0]
            publication_date = datetime.date(int(float(tmp['Year'])), int(float(tmp['Month'])), int(float(tmp['Day'])))
        else:
            publication_date = None

        corresponding_authors.append({
            'first_name': author['ForeName'] if 'ForeName' in author else None,
            'last_name': author['LastName'] if 'LastName' in author else None,
            'email': email,
            'affiliation': affiliation,
            'pubmed_id': pubmed_id,
            'n_citation': n_citation,
            'abstract': abstract_text,
            'publication_date': publication_date,
            })

    return corresponding_authors

if __name__ == '__main__':
    main()
