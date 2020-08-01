# geneutils

[![PyPI version](https://badge.fury.io/py/geneutils.svg)](https://badge.fury.io/py/geneutils)

This tool was originally inspired by the [work here](https://github.com/Gurdhhu/bioinf_scripts) for dealing with nucleotide sequences and more tools will be added as needed. This toolbox requires python 3 and can be installed.

## Table of contents
* [Installation](#installation)
* [geneutils cli tools](#geneutils-cli-tools)
    * [blasthit](#blasthit)


## Installation
This assumes that you have native python 3 & pip installed in your system, you can test this by going to the terminal (or windows command prompt) and trying

```python``` and then ```pip list```

If you get no errors and you have python 3 or higher you should be good to go. To install **geneutils** you can install using two methods

```pip install geneutils```

or you can also try

```pip install geneutils --user```

or you can also try

```
git clone https://github.com/samapriya/geneutils.git
cd geneutils
python setup.py install
```

Installation is an optional step; the application can be also run directly by executing pydrop.py script. The advantage of having it installed is being able to execute porg as any command line tool. I recommend installation within virtual environment. If you don't want to install, browse into the geneutils folder and try ```python geneutils.py -h``` to get to the same result.

## geneutils cli tools
This is a command line tool and it is designed to simply call the tools you need.

![geneutils_main](https://user-images.githubusercontent.com/6677629/89102487-b2d5b200-d3d7-11ea-937e-cd6e661de31c.gif)

### blasthit
This script is intended for taxonomic annotation of **blast** results (blastn, tblastn, blastp or blastx) saved in **Hit Table CSV** format where GenBank accession numbers are in the 4th column. It uses **efetch** function from **Bio.Entrez** package to get information about accessions from **GenBank Nucleotide (Nuccore)** or **Protein** databases.
The output is an annotated CSV file "*_annotated.csv" with the following columns added:

* Record name
* Species name
* Full taxonomy
* Reference
* Date of update

![geneutils_bhits](https://user-images.githubusercontent.com/6677629/89102488-b701cf80-d3d7-11ea-8f1b-b26886563e84.gif)

| arguments | description |
| --- | --- |
| path | Pathway to csv-formatted Hit Table file with blastn results. Positional argument |
| db | For the output of **nucleotide blast** or **tblastn**, use <code>n</code>. For the output of **protein blast** or **blastx**, use <code>p</code>. Positional argument |
| email | Your NCBI email. Positional argument |
| -h, --help | Show help message and exit. Optional argument |
