""" Fetch author emails from Pubmed

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-12
:Copyright: 2016, Karr Lab
:License: MIT
"""

from Bio import Entrez
import csv

# configure entrez
Entrez.email = 'karr@mssm.edu'
query = 'systems biology'
max_results = 10000


def parse_article_metadata(article):
    medline_article = article['MedlineCitation']['Article']
    if 'AuthorList' not in medline_article:
        return None

    author_list = medline_article['AuthorList']
    if not author_list:
        return None

    author = author_list[0]
    affiliation_info = author['AffiliationInfo']
    if not affiliation_info:
        return None

    affiliation = affiliation_info[0]['Affiliation']
    affiliation, _, email = affiliation[0:-1].rpartition(' ')
    if '@' not in email:
        return None

    return (author['LastName'] if 'LastName' in author else None,
            author['ForeName'] if 'ForeName' in author else None,
            email,
            affiliation,
            )


# search PubMed
search_handle = Entrez.esearch(db='pubmed',
                               sort='relevance',
                               retmax=str(max_results),
                               term=query)
search_results = Entrez.read(search_handle)
pubmed_ids = search_results['IdList']

# fetch publication metadata
people = set()
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
        person = parse_article_metadata(article)
        if person:
            people.add(person)

# save author information to a tsv file
with open('emails.tsv', 'w') as file:
    csv_writer = csv.writer(file, delimiter='\t')
    csv_writer.writerow(('First', 'Last', 'Email', 'Affiliation'))
    for person in people:
        csv_writer.writerow(person)
