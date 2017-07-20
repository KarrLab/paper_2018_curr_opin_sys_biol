#!/usr/bin/env python

""" Fetch author emails from Europe PMC

See also the `documentation for the Europe PMC webservice <https://europepmc.org/RestfulWebService#search>`_

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-07-19
:Copyright: 2016, Karr Lab
:License: MIT
"""

from six import StringIO
import codecs
import csv
import datetime
import json
import math
import requests
import six

EUROPE_PMC_ENDPOINT = 'http://www.ebi.ac.uk/europepmc/webservices/rest/search'
PAGE_SIZE = 10  # maximum is 1000
MAX_PAGES = 2#float('inf')

query = (
    '('
    '    (TITLE:"systems biology" OR ABSTRACT:"systems biology" OR KW:"systems biology") OR '
    '    (TITLE:"computational neuroscience" OR ABSTRACT:"computational neuroscience" OR KW:"Models, Neurological") OR '
    '    (TITLE:"multiscale modeling" OR ABSTRACT:"multiscale modeling")'
    ') AND ('
    '          ('
    '              (TITLE:model OR ABSTRACT:model OR KW:"Models, Biological") AND '
    '              (TITLE:MATHEMATICAL OR ABSTRACT:mathematical OR ABSTRACT:computational OR TITLE:computational)'
    '          ) OR '
    '    (TITLE:simulation OR ABSTRACT:simulation OR KW:"Computer Simulation")'
    ')AND ('
    '    FIRST_PDATE:[2012-01-01 TO 2020-12-31]'
    ')'
)


class HashableDict(dict):

    def __hash__(self):
        return hash(tuple(sorted(self.items())))


class UnicodeWriter:
    """
    A CSV writer which will write rows to CSV file "f",
    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        # Redirect output to a queue
        self.queue = StringIO()
        self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
        self.stream = f
        self.encoder = codecs.getincrementalencoder(encoding)()

    def writerow(self, row):
        self.writer.writerow([s.encode("utf-8") for s in row])
        # Fetch UTF-8 output from the queue ...
        data = self.queue.getvalue()
        data = data.decode("utf-8")
        # ... and reencode it into the target encoding
        data = self.encoder.encode(data)
        # write to the target stream
        self.stream.write(data)
        # empty queue
        self.queue.truncate(0)

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""

    if obj is None:
        return ''
    if isinstance(obj, datetime.date):
        return str(obj)
    raise TypeError("Type %s not serializable" % type(obj))


def get_email(affiliation):
    if '@' not in affiliation:
        return None

    for word in reversed(affiliation.split(' ')):
        if '@' in word:
            if word[-1] == '.':
                return word[0:-1]
            else:
                return word


def merge_investigators(a, b):
    a['id'].update(b['id'])
    a['email'].update(b['email'])
    a['articles'].update(b['articles'])
    return a


article_citation_counts = []
nextCursorMark = '*'
page = 0

investigators_by_email = {}
investigators_by_id = {}
while True:
    response = requests.get(EUROPE_PMC_ENDPOINT, params={
        'query': query,
        'resulttype': 'core',
        'format': 'json',
        'pageSize': PAGE_SIZE,
        'cursorMark': nextCursorMark,
    })
    response.raise_for_status()
    response_data = response.json()

    page += 1
    n_pages = int(math.ceil(float(response_data['hitCount']) / float(PAGE_SIZE)))
    print('Downloading page {} of {} ...'.format(page, n_pages))

    for article in response_data['resultList']['result']:
        if 'authorList' not in article:
            continue

        if 'author' not in article['authorList']:
            continue

        for author in article['authorList']['author']:
            if 'authorId' in author:
                id = HashableDict(type=author['authorId']['type'], value=author['authorId']['value'])
            else:
                id = None

            if 'affiliation' in author:
                affiliation = author['affiliation']
                email = get_email(author['affiliation'])
            else:
                affiliation = None
                email = None

            if not id and not email:
                continue

            first_name = author['firstName'] if 'firstName' in author else None
            last_name = author['lastName'] if 'lastName' in author else None

            title = article['title']
            if 'journalInfo' in article and 'journal' in article['journalInfo'] and 'title' in article['journalInfo']['journal']:
                journal = article['journalInfo']['journal']['title']
            else:
                journal = ''
            year, month, day = article['firstPublicationDate'].split('-')
            date = datetime.date(int(float(year)), int(float(month)), int(float(day)))
            doi = article['doi'] if 'doi' in article else None
            pmid = article['pmid'] if 'pmid' in article else None
            pmcid = article['pmcid'] if 'pmcid' in article else None

            citations = article['citedByCount'] if 'citedByCount' in article else math.nan

            investigator = {
                'id': set([id]) if id else set([]),
                'email': set([email]) if email else set([]),
                'articles': set([HashableDict(
                    first_name=first_name,
                    last_name=last_name,
                    email=email,
                    id=id,
                    affiliation=affiliation,
                    title=title,
                    journal=journal,
                    date=date,
                    doi=doi,
                    pmid=pmid,
                    pmcid=pmcid,
                    citations=citations,
                )]),
            }

            if email and email in investigators_by_email:
                investigator = merge_investigators(investigators_by_email[email], investigator)
            if id and id in investigators_by_id:
                investigator = merge_investigators(investigators_by_id[id], investigator)

            for email in investigator['email']:
                investigators_by_email[email] = investigator
            for id in investigator['id']:
                investigators_by_id[id] = investigator

    if page == n_pages or page >= MAX_PAGES:
        break

    nextCursorMark = response_data['nextCursorMark']

# get list of investigators
investigators = []
for investigator in investigators_by_email.values():
    if investigator not in investigators:
        investigators.append(investigator)
for investigator in investigators_by_id.values():
    if investigator not in investigators:
        investigators.append(investigator)


def key_func(investigator):
    latest_article = sorted(investigator['articles'], key=lambda article: article['date'], reverse=True)[0]
    return (latest_article['last_name'], latest_article['first_name'])
investigators.sort(key=key_func)

# save author information to a tsv file
with open('investigators.tsv', 'w') as file:
    if six.PY2:
        csv_writer = UnicodeWriter(file, delimiter='\t')
    else:
        csv_writer = csv.writer(file, delimiter='\t')
    csv_writer.writerow([
        'Last name', 'First name', 'Email', 'All emails', 'All IDs',
        'Affiliation', 'All affiliations',
        'Lastest article title', 'Lastest article journal', 'Lastest article date',
        'Lastest article DOI', 'Lastest article PMID', 'Lastest article PMCID',
        'Number articles', 'Articles', 'Citations',
    ])
    for investigator in investigators:
        if not investigator['email']:
            continue
        latest_article = sorted(investigator['articles'], key=lambda article: article['date'], reverse=True)[0]

        latest_email = sorted(filter(lambda article: article['email'] is not None, investigator['articles']), key=lambda article: article['date'], reverse=True)[0]['email']

        csv_writer.writerow([
            latest_article['last_name'] or '',
            latest_article['first_name'] or '',
            latest_email or '',
            json.dumps(list(investigator['email'])) or '',
            json.dumps(list(investigator['id'])) or '',
            latest_article['affiliation'] or '',
            json.dumps(list(set([article['affiliation'] for article in investigator['articles']]))) or '',
            latest_article['title'] or '',
            latest_article['journal'] or '',
            str(latest_article['date']) or '',
            latest_article['doi'] or '',
            latest_article['pmid'] or '',
            latest_article['pmcid'] or '',
            str(len(investigator['articles'])),
            json.dumps(list(investigator['articles']), default=json_serial) or '',
            str(sum(filter(lambda citations: not math.isnan(citations), [article['citations'] for article in investigator['articles']]))),
        ])

# print status message
print('Done!')
