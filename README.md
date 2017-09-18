# Code for the Karr Lab's 2018 review in *Current Opinion in Systems Biology*
This repository provides the code for the Karr Lab's 2018 whole-cell modeling review in *Current Opinion in Systems Biology*.

The code performs the following functions:

* Generates part of the email list used for the biomodeling community survey reported in the paper
  * Uses the `Europe PMC webservice <https://europepmc.org/RestfulWebService#search>`_ to find articles about biomodeling
  * Extracts the author emails and ORCIDs and citation counts from the returned results
  * Groups the articles by author (common email and/or ORCID)
  * Counts the number of articles and sums the citations for each author
  * Saves the results to a TSV file
* Generates the supplementary figures

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

The code will save the results to `investigators.tsv` and `figures/`.

## How to cite this code
Check back in 2018

## License
The code is released under the [MIT license](LICENSE).

## Development team
This code was developed by the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.
