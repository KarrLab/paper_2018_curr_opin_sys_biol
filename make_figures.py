""" Plot survey results for supplementary material

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-09-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

from matplotlib import pyplot
from matplotlib.backends import backend_pdf
import itertools
import numpy
import os
import pyexcel_xlsx

invitations_filename = os.path.join(os.path.dirname(__file__), 'survey_invitations.xlsx')
# path to Excel workbook with survey invitation list

responses_filename = os.path.join(os.path.dirname(__file__), 'survey_responses.xlsx')
# path to Excel workbook with survey responses

OUT_DIR = os.path.join(os.path.dirname(__file__), 'figures')
# output directory

#################################
# make the output directory if it doesn't already exist
#################################
if not os.path.isdir(OUT_DIR):
    os.makedirs(OUT_DIR)
    
#################################
# read invitations list
#################################
rows = pyexcel_xlsx.get_data(invitations_filename)['MASTER_emails']
invitations = []
for row in rows[1:]:
    invitation = {}
    invitation['email'] = row[0]
    invitation['first'] = row[1]
    invitation['last'] = row[2]
    invitation['affilitation'] = row[3]
    invitation['position'] = row[4]
    invitation['sector'] = row[5]
    invitation['country'] = row[6]
    invitation['invitation_date'] = row[7]
    
    if len(row) >= 9:    
        invitation['bounced'] = row[8] == 'Y'
    else:
        invitation['bounced'] = False
    
    if len(row) >= 10:
        invitation['responded'] = row[9] == 'Y'
    else:
        invitation['responded'] = False
    
    invitations.append(invitation)

#################################
# read survey responses
#################################
# todo
rows = pyexcel_xlsx.get_data(responses_filename)['Survey responses']
responses = []
for row in rows[1:]:
    response = {
        'exclude': row[25],
        'fields': row[26].split(';'),
        'subfield': row[27].split(';'),
        'bio_focus': row[28].split(';'),
        'model_size': row[29].split(';'),
        'data_types': row[30].split(';'),
        'data_sources': row[31].split(';'),
        'formalisms': row[32].split(';'),
        'tools': row[33].split(';'),
        'model_formats': row[34].split(';'),
        'repositories': row[35].split(';'),
        'programming_languages': row[36].split(';'),
        'other_tools': row[37].split(';'),
        'bottlenecks': row[38].split(';'),
        'key_problems': row[39].split(';'),
        'key_bottlenecks': row[40].split(';'),
        'needed_resources': row[41].split(';'),
        'needed_meetings': row[42].split(';'),
        'other_thoughts': row[43].split(';'),
        'sector': row[44],
        'country': row[45],
    }
    responses.append(response)
    
responses = list(filter(lambda r: r['exclude'] != 'Y', responses))

#################################
# make and save figures
#################################    
def plot_histogram(groups, x_label, y_label, filename, sort=True):
    if sort:
        groups.sort(key=lambda group: -1 if group['label'] == 'Unknown' else group['count'])
    labels = [group['label'] for group in groups]
    counts = [group['count'] for group in groups]
    
    # plot data
    fig, axes = pyplot.subplots(nrows=1, ncols=1)    
    y = numpy.arange(len(groups))
    axes.barh(y, counts)
    axes.set_xlim([0, numpy.max(counts)])
    axes.set_ylim([-0.5, len(groups)])
    axes.set_yticks(y + 0.5)
    axes.set_yticklabels(labels)    
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    
    #todo: label percent
    #todo: change color
    #todo: set bar size, page size, font size
    
    # save figure
    with backend_pdf.PdfPages(os.path.join(OUT_DIR, filename)) as pp:
        pp.savefig(fig, transparent=True, bbox_inches='tight', pad_inches=0.)
    pyplot.close(fig)
    
def plot_pie(groups, filename):
    groups.sort(key=lambda group: -1 if group['label'] == 'Unknown' else group['count'])
    labels = [group['label'] for group in groups]
    counts = [group['count'] for group in groups]
    
    # plot data
    fig, axes = pyplot.subplots(nrows=1, ncols=1)        
    explode = numpy.full((len(counts),), 0.01)
    axes.pie(counts, explode=explode, labels=labels, autopct='%d%%', shadow=False, startangle=90)
    axes.axis('equal')
    
    # save figure
    with backend_pdf.PdfPages(os.path.join(OUT_DIR, filename)) as pp:
        pp.savefig(fig, transparent=True, bbox_inches='tight', pad_inches=0.)
    pyplot.close(fig)

# invitations by country
groups = []
key = lambda i: i['country'] or 'Unknown'
invitations.sort(key=key)
for key, group in itertools.groupby(invitations, key):
    groups.append({'label': key, 'count': len(list(group))})
plot_histogram(groups, 'Number of researchers', 'Country', 'Invitations by country.pdf')

# invitations by response
groups = [
    {'label': 'Yes', 'count': 0}, 
    {'label': 'No', 'count': 0},
    ]
for invitation in invitations:
    if invitation['responded']:
        groups[0]['count'] += 1
    elif not invitation['bounced']:
        groups[1]['count'] += 1
plot_histogram(groups, 'Number of researchers', 'Responded to survey', 'Invitations by response.pdf', sort=False)

# responses by sector
groups = []
key = lambda r: r['sector'] or 'Unknown'
responses.sort(key=key)
for key, group in itertools.groupby(invitation_responses, key):
    groups.append({'label': key, 'count': len(list(group))})
plot_histogram(groups, 'Number of researchers', 'Sector', 'Responses by sector.pdf')

# responses by country
groups = []
key = lambda r: r['country'] or 'Unknown'
responses.sort(key=key)
for key, group in itertools.groupby(responses, key):
    groups.append({'label': key, 'count': len(list(group))})
plot_histogram(groups, 'Number of researchers', 'Country', 'Responses by country.pdf')

