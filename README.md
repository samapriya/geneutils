# geneutils

[![PyPI version](https://badge.fury.io/py/geneutils.svg)](https://badge.fury.io/py/geneutils)
![CI geneutils](https://github.com/samapriya/geneutils/workflows/CI%20geneutils/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3975287.svg)](https://doi.org/10.5281/zenodo.3975287)

This tool was originally inspired by the [work here](https://github.com/Gurdhhu/bioinf_scripts) for dealing with nucleotide sequences and more tools will be added as needed. This toolbox requires python 3 and can be installed.

## Table of contents
* [Installation](#installation)
* [geneutils cli tools](#geneutils-cli-tools)
    * [init](#init)
    * [blasthit](#blasthit)


## Installation
This assumes that you have native python 3 & pip installed in your system, you can test this by going to the terminal (or windows command prompt) and trying

```python``` and then ```pip list```

If you get no errors and you have python 3 or higher you should be good to go. To install **geneutils** you can install using two methods

```pip install geneutils```

or you can also try

```pip install geneutils --user```


To upgrade to a new version from old you can try

```
pip install geneutils --upgrade
```

or

```
pip install geneutils --user --upgrade
```

or you can also try

```
git clone https://github.com/samapriya/geneutils.git
cd geneutils
python setup.py install
```

Though there is a conda installer for linux/mac for geneutils

```conda install -c samapriya geneutils```

The recommended way would be to use the conda environment/terminal and then do a ```pip install geneutils```

Installation is an optional step; the application can be also run directly by executing geneutils.py script. The advantage of having it installed is being able to execute porg as any command line tool. I recommend installation within virtual environment. If you don't want to install, browse into the geneutils folder and try ```python geneutils.py -h``` to get to the same result.

## geneutils cli tools
This is a command line tool and it is designed to simply call the tools you need.

![geneutils_main](https://user-images.githubusercontent.com/6677629/89606925-3f9cc780-d83f-11ea-9a05-6b1c9a31ff68.gif)

### init
Turns out there are benefits of registering for a NCBI account and to use your email address and your API key. Apart from NCBI having a way of contacting you, the API key raises the rate limit imposed on your queries. This is a recommended step and a user should not skip thought it is possible to use the blashit tool without the API key.
*From the NCBI account page you can find this information*

**E-utils users are allowed 3 requests/second without an API key. Create an API key to increase your e-utils limit to 10 requests/second......Only one API Key per user. Replacing or deleting will inactivate the current key. Use this key by passing it with api_key=API_KEY parameter.**

Generate your API Key
![ncbi_apikey](https://user-images.githubusercontent.com/6677629/89606628-67d7f680-d83e-11ea-9c43-328903dcd6b7.gif)

The init tool saves your email address and API key to be saved in your local machine which can be used instead of typing out your email over and over again. The API key is a clear entry meaning you cannot see when you type in or paste your API key for safety.
![ncbi_cred](https://user-images.githubusercontent.com/6677629/89147373-6c419e00-d524-11ea-8043-58f3e9699b5f.gif)

### blasthit
This script is intended for taxonomic annotation of **blast** results (blastn, tblastn, blastp or blastx) saved in **Hit Table CSV** format where GenBank accession numbers are in the 4th column. It uses **efetch** function from **Bio.Entrez** package to get information about accessions from **GenBank Nucleotide (Nuccore)** or **Protein** databases.
The output is an annotated CSV file "*_annotated.csv" with the following columns added:

* Record name
* Species name
* Full taxonomy
* Reference
* Date of update

![geneutils_bhits](https://user-images.githubusercontent.com/6677629/89607253-0dd83080-d840-11ea-997b-9b69cbb4e8b4.gif)

| arguments | description |
| --- | --- |
| path | Pathway to csv-formatted Hit Table file with blastn results. Positional argument |
| db | For the output of **nucleotide blast** or **tblastn**, use <code>n</code>. For the output of **protein blast** or **blastx**, use <code>p</code>. Positional argument |
| email | Your NCBI email. Optional argument |
| -h, --help | Show help message and exit. Optional argument |


## Changelog

### 0.0.4
  - Fixed issue with path basename
  
### 0.0.3
  - Fixed issue with repeating email id for accession blocks
  - Updated ReadMe to include **geneutils init** as first step.
  - Fixed CSV write issue outside loop

### 0.0.2
  - Added credential tool to save NCBI email and API Key
  - Minor fixes to overall functionality.
