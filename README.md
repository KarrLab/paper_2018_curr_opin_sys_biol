# Code for "A blueprint for human whole-cell modeling"

This repository provides the code for "A blueprint for human whole-cell modeling" by Szigeti et al.

The code performs the following functions:

* Generates part of the email list used for the biomodeling community survey reported in the paper
  * Uses the [Europe PMC webservice](https://europepmc.org/RestfulWebService#search) to find articles about biomodeling
  * Extracts the author emails and ORCIDs and citation counts from the returned results
  * Groups the articles by author (common email and/or ORCID)
  * Counts the number of articles and sums the citations for each author
  * Saves the results to a TSV file
* Analyzes the invited investigators and the survey responses
  * Generates several supplementary figures

## How to install this code
1. Install Python and pip
    ```
    apt-get install python python-pip
    ````
3. Install the require PyPI packages
    ```
    pip install -r requirements.txt
    ````
  
## How to run this code
```
python get_author_emails.py
python make_figures.py
```

The first script generates part of the list of investigators that we invited to participate in the survey and saves this list to `investigators.tsv`. The second script analyzes the survey invitations in `survey_invitations.xlsx` and survey responses in `survey_responses-edited.xlsx` and saves several figures to `figures/`.

## How to cite this code
Szigeti B, Roth YD, Sekar JAP, Goldberg AP, Pochiraju SC & Karr JR†. A blueprint for human whole-cell modeling. Curr Opin Syst Biol 7, 8–15 (2018). doi: [10.1016/j.coisb.2017.10.005](https://doi.org/10.1016/j.coisb.2017.10.005)

## License
The code is released under the [MIT license](LICENSE).

## Development team
This code was developed by the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.
