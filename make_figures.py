""" Plot survey results for supplementary material

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2017-09-08
:Copyright: 2017, Karr Lab
:License: MIT
"""

from matplotlib import pyplot
from matplotlib.backends import backend_pdf
import itertools
import matplotlib
import natsort
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
rows = pyexcel_xlsx.get_data(invitations_filename)['Invitations']
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


def parse_list(row, col):
    if len(row) > col:
        str = row[col]
    else:
        return []

    if str:
        return str.split(';')
    return []

rows = pyexcel_xlsx.get_data(responses_filename)['Responses']
responses = []
for row in rows[2:]:
    if row[25] == 'Y':
        continue

    response = {
        'sectors': parse_list(row, 26),
        'fields': parse_list(row, 27),
        'countries': parse_list(row, 28),
        'focus': parse_list(row, 29),
        'bio_focus': parse_list(row, 30),
        'model_size': parse_list(row, 31),
        'data_types': parse_list(row, 32),
        'data_sources': parse_list(row, 33),
        'formalisms': parse_list(row, 34),
        'tools': parse_list(row, 35),
        'model_formats': parse_list(row, 36),
        'repositories': parse_list(row, 37),
        'programming_languages': parse_list(row, 38),
        'other_tools': parse_list(row, 39),
        'bottlenecks': parse_list(row, 40),
        'key_bottlenecks': parse_list(row, 41),
        'key_problems': parse_list(row, 42),
        'needed_resources': parse_list(row, 43),
        'needed_meetings': parse_list(row, 44),
        'other_thoughts': parse_list(row, 45),
    }
    responses.append(response)

#################################
# make and save figures
#################################


def group_invitations(invitations, key_name):
    key = lambda i: i[key_name] or 'Unknown'
    invitations.sort(key=natsort.natsort_keygen(key, alg=natsort.IGNORECASE))

    groups = []
    for key, group in itertools.groupby(invitations, key):
        groups.append({'label': key, 'count': len(list(group))})

    return groups


def group_responses(responses, key_name, unknown=True):
    all_responses = []
    for response in responses:
        if response[key_name]:
            all_responses += response[key_name]
        elif unknown:
            all_responses.append('Unknown')

    key = lambda el: el
    all_responses.sort(key=natsort.natsort_keygen(key, alg=natsort.IGNORECASE))
    groups = []
    for key, group in itertools.groupby(all_responses, key):
        groups.append({'label': key, 'count': len(list(group))})

    return groups


def plot_histogram(groups, total, x_label, y_label, filename, threshold_count=1.):
    def sort_func(group):
        if group['label'] == 'Other':
            return (1, group['label'])
        elif group['label'] == 'Unknown':
            return (2, group['label'])
        else:   
            return (-1 * float(group['count']), group['label'])

    groups.sort(key=natsort.natsort_keygen(sort_func, alg=natsort.IGNORECASE), reverse=True)
    labels = [group['label'] for group in groups]
    counts = numpy.array([group['count'] for group in groups])

    filtered_labels = []
    filtered_counts = []
    other_count = 0.
    other_caption = []
    for label, count in zip(labels, counts):
        if count >= threshold_count:
            filtered_labels.append(label)
            filtered_counts.append(count)
        else:
            other_count += count
            other_caption.append('{} ({})'.format(label, count))
    labels = filtered_labels
    counts = filtered_counts
    if other_count:
        if len(labels) > 0 and labels[0] == 'Unknown':
            labels.insert(1, 'Other')
            counts.insert(1, other_count)
        else:
            labels.insert(0, 'Other')
            counts.insert(0, other_count)
        other_caption.reverse()
        if len(other_caption) == 1:        
            caption = 'Other: ' + other_caption[-1] + '.'
        else:
            caption = 'Other: ' + ', '.join(other_caption[0:-1]) + ' and ' + other_caption[-1] + '.'
    else:
        caption = ''

    # plot data    
    h_bar = 0.12
    h_axis = 0.35
    h_fig = (len(labels) + 0.1) * h_bar + h_axis

    w_fig = 6.5
    w_axis = 1.68
    w_label = 0.30
    w_grid = w_fig - w_axis - w_label

    fig = pyplot.figure()
    fig.set_size_inches(w_fig, h_fig)
    fig.set_dpi(300.)
    fig.set_frameon(False)
    axes = fig.add_axes((w_axis / w_fig, h_axis / h_fig, w_grid / w_fig, 1 - h_axis / h_fig))

    y = numpy.arange(len(labels))
    axes.barh(y, counts, left=0, tick_label=labels, height=0.9, color='black')
    axes.set_xlim([0, numpy.max(counts)])
    axes.set_ylim([-0.55, len(labels) - 0.45])
    for tick in axes.xaxis.get_major_ticks():
        tick.label.set_fontsize(7)
        tick.label.set_fontname('Arial')
    for tick in axes.yaxis.get_major_ticks():
        tick.label.set_fontsize(7)
        tick.label.set_fontname('Arial')
    axes.set_xlabel(x_label, fontsize=9, fontname='Arial')
    axes.set_ylabel(y_label, fontsize=9, fontname='Arial')

    for label, count, rect in zip(labels, counts, axes.patches):
        width = rect.get_width()
        axes.text(rect.get_x() + rect.get_width() + 0.02 / w_grid * numpy.max(counts), rect.get_y() + rect.get_height() / 2, 
            '{:.1f}%'.format(float(count) / float(total) * 100.), ha='left', va='center', fontsize=7, fontname='Arial')

    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)

    # save figure
    with backend_pdf.PdfPages(os.path.join(OUT_DIR, filename + '.pdf')) as pp:
        pp.savefig(fig, transparent=True)

    pyplot.close(fig)

    return caption

figures = [
    {
        'data': 'invitations',
        'title': 'Sectors of the invited scientists.',
        'short_title': 'Sectors of the invited scientists.',
        'key': 'sector',
        'y_axis_label': 'Sector',
        'threshold_count': 0.,
        'unknown': True,
    },
    {
        'data': 'invitations',
        'title': 'Countries of the invited scientists.',
        'short_title': 'Countries of the invited scientists.',
        'key': 'country',
        'y_axis_label': 'Country',
        'threshold_count': 5.,
        'unknown': True,
    },
    {
        'data': 'responses',
        'title': 'What sector do you work in?',
        'short_title': 'What sector do you work in?',
        'key': 'sectors',
        'y_axis_label': 'Sector',
        'threshold_count': 0.,
        'unknown': True,
    },
    {
        'data': 'responses',
        'title': 'What is your research field?',
        'short_title': 'What is your research field?',
        'key': 'fields',
        'y_axis_label': 'Field',
        'unknown': True,
        'threshold_count': 3.,
    },
    {
        'data': 'responses',
        'title': 'Where do you work?',
        'short_title': 'Where do you work?',
        'key': 'countries',
        'y_axis_label': 'Country',
        'threshold_count': 3.,
        'unknown': True,
    },
    {
        'data': 'responses',
        'title': 'What is the primary focus on your research?',
        'short_title': 'What is the primary focus on your research?',
        'key': 'focus',
        'y_axis_label': 'Research focus',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What is the biological focus of your research?',
        'short_title': 'What is the biological focus of your research?',
        'key': 'bio_focus',
        'y_axis_label': 'Biological focus',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, how large are the models that you typically use in your research?',
        'short_title': 'How large are the models that you typically use in your research?',
        'key': 'model_size',
        'y_axis_label': 'Model size',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, what types of data sources do you use?',
        'short_title': 'What types of data sources do you use?',
        'key': 'data_types',
        'y_axis_label': 'Data type',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, which data sources do you use?',
        'short_title': 'Which data sources do you use?',
        'key': 'data_sources',
        'y_axis_label': 'Data source',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, which mathematical representations and simulation algorithms do you use most frequently?',
        'short_title': 'Which mathematical representations do you use most frequently?',
        'key': 'formalisms',
        'y_axis_label': 'Mathematical formalism',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, which tools do you most frequently use to build and/or simulate models?',
        'short_title': 'Which tools do you most frequently use to build and/or simulate models?',
        'key': 'tools',
        'y_axis_label': 'Tool',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, which languages do you most frequently use to represent models?',
        'short_title': 'Which languages do you most frequently use to represent models?',
        'key': 'model_formats',
        'y_axis_label': 'Model format',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'If you use models in your research, which resources do you most frequently use to distribute models?',
        'short_title': 'Which resources do you most frequently use to distribute models?',
        'key': 'repositories',
        'y_axis_label': 'Resource',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'Which programming languages do you most frequently use in your research?',
        'short_title': 'Which programming languages do you most frequently use?',
        'key': 'programming_languages',
        'y_axis_label': 'Programming language',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What additional tools do you frequently use in your research?',
        'short_title': 'What additional tools do you frequently use in your research?',
        'key': 'other_tools',
        'y_axis_label': 'Tool',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What are the most time-consuming aspects of your research?',
        'short_title': 'What are the most time-consuming aspects of your research?',
        'key': 'bottlenecks',
        'y_axis_label': 'Bottleneck',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What do think are the main bottlenecks to building more predictive models?',
        'short_title': 'What do think are the main bottlenecks to building more predictive models?',
        'key': 'key_bottlenecks',
        'y_axis_label': 'Bottleneck',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'Given access to improved data and tools, what biomedical problems would you like to use models to solve (e.g. designer microorganisms, personalized cancer therapy)?',
        'short_title': 'What biomedical problems would you like to use models to solve?',
        'key': 'key_problems',
        'y_axis_label': 'Problem',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What additional databases, methods, tools, and/or standards would help accelerate your research?',
        'short_title': 'What additional resources would help accelerate your research?',
        'key': 'needed_resources',
        'y_axis_label': 'Resource',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'What additional meetings, courses, and/or training opportunities would help accelerate your research?',
        'short_title': 'What additional meetings would help accelerate your research?',
        'key': 'needed_meetings',
        'y_axis_label': 'Meeting',
        'threshold_count': 3.,
        'unknown': False,
    },
    {
        'data': 'responses',
        'title': 'Are there any additional thoughts that you would like to share?',
        'short_title': 'Are there any additional thoughts that you would like to share?',
        'key': 'other_thoughts',
        'y_axis_label': 'Other thought',         
        'threshold_count': 3.,
        'unknown': False,
    },
]

for i_figure, figure in enumerate(figures):
    if figure['data'] == 'invitations':
        groups = group_invitations(invitations, figure['key'])
        total = len(invitations)
        x_axis_label = 'Number of invited scientists'
        filename = '{} - Invitations by {}'.format(i_figure + 1, figure['y_axis_label'].lower())

    else:
        groups = group_responses(responses, figure['key'], unknown=figure['unknown'])
        total = len(responses)
        x_axis_label = 'Number of responses'
        filename = '{} - Responses by {}'.format(i_figure + 1, figure['y_axis_label'].lower())

    figure['filename'] = filename
    figure['caption'] = plot_histogram(groups, total, x_axis_label,
                                       figure['y_axis_label'], filename,
                                       threshold_count=figure['threshold_count'])

# print figures
with open(os.path.join(OUT_DIR, 'figures.tex'), 'w') as file:
    for figure in figures:
        file.write('\clearpage\n\suppfig{{{}}}{{{}}}{{{}}}{{{}}}\n\n'.format(
            figure['filename'],
            figure['title'],
            figure['short_title'],
            figure['caption'].replace('&', '\&').replace('#', '\#').replace('%', '\%'),
        ))
