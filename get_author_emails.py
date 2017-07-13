#!/usr/bin/env python

""" Fetch author emails from Pubmed

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-12
:Copyright: 2016, Karr Lab
:License: MIT
"""

from Bio import Entrez
import csv
import datetime
import itertools
import math
import numpy

# configure entrez
Entrez.email = 'karr@mssm.edu'
query = (
    '(simulation OR ((mathematical OR computational) AND model))'
    ' AND ("systems biology" OR multi-scale OR multiscale OR "computational neuroscience")'
    ' AND ("2012/01/01"[Date - Publication] : "2099/12/31"[Date - Publication])'
    )
max_results = 10000


def parse_article_metadata(article):
    medline_article = article['MedlineCitation']['Article']
    if 'AuthorList' not in medline_article:
        return []

    author_list = medline_article['AuthorList']
    if not author_list:
        return []

    corresponding_authors = []
    for author in author_list:
        affiliation_info = author['AffiliationInfo']
        if not affiliation_info:
            continue

        affiliation = affiliation_info[0]['Affiliation']
        affiliation, _, email = affiliation[0:-1].rpartition(' ')
        if '@' not in email:
            continue

        id_list = article['PubmedData']['ArticleIdList']
        pubmed_id = None
        for id in id_list:
            if id.attributes['IdType'] == 'pubmed':
                pubmed_id = str(id)

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
            'publication_date': publication_date,
            })

    return corresponding_authors


# search PubMed
search_handle = Entrez.esearch(db='pubmed',
                               sort='relevance',
                               retmax=str(max_results),
                               term=query)
search_results = Entrez.read(search_handle)
pubmed_ids = search_results['IdList']

# fetch publication metadata
corresponding_authors = []
chunk_size = 1000
for i_chunk in range(int(math.ceil(len(pubmed_ids) / chunk_size))):
    print('Retrieving articles {}..{} of {}'.format(
        i_chunk * chunk_size + 1,
        min(len(pubmed_ids), (i_chunk + 1) * chunk_size),
        len(pubmed_ids)))

    fetch_handle = Entrez.efetch(db='pubmed', retmode='xml', id=','.join(
        pubmed_ids[i_chunk * chunk_size:min(len(pubmed_ids), (i_chunk + 1) * chunk_size)]))
    article_metadata = Entrez.read(fetch_handle)

    # extract author info from results
    for article in article_metadata['PubmedArticle']:
        corresponding_authors.extend(parse_article_metadata(article))

# group authors by email
investigators = []
key = lambda author: author['email'].lower()
corresponding_authors.sort(key=key)
for email, articles in itertools.groupby(corresponding_authors, key):
    sorted_articles = sorted(articles, 
        key=lambda article: article['publication_date'] if article['publication_date'] else datetime.date(1, 1, 1), 
        reverse=True)
    lastest_article = sorted_articles[0]

    investigators.append({
        'first_name': lastest_article['first_name'],
        'last_name': lastest_article['last_name'],
        'email': email,
        'affiliation': lastest_article['affiliation'],
        'all_affiliations': [article['affiliation'] for article in sorted_articles],
        'lastest_publication_date': lastest_article['publication_date'],
        'pubmed_ids': [article['pubmed_id'] for article in sorted_articles],
        })

# sort investigators by name
investigators.sort(key=lambda investigator: (investigator['last_name'], investigator['first_name']))

# count number of authors 
article_count_frequencies = numpy.bincount([len(i['pubmed_ids']) for i in investigators])

print('{:12}  {:17}'.format('No. articles', 'No. investigators'))
for n_articles, n_investigators in enumerate(article_count_frequencies):
    print('{:12}  {:17}'.format(n_articles, n_investigators))

# save author information to a tsv file
with open('investigators.tsv', 'w') as file:
    csv_writer = csv.DictWriter(file, delimiter='\t', fieldnames=investigators[0].keys())
    csv_writer.writeheader()
    for investigator in investigators:
        csv_writer.writerow(investigator)
