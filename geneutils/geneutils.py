#!/usr/bin/python
# -*- coding: utf-8 -*-

__copyright__ = """

MIT License

Copyright (c) 2020 Samapriya Roy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
__license__ = "MIT License"

import csv
import sys
import os
import json
import getpass
import requests
import pkg_resources
import pandas as pd
import argparse
from bs4 import BeautifulSoup
from os.path import expanduser
from Bio import Entrez


class Solution:
    def compareVersion(self, version1, version2):
        versions1 = [int(v) for v in version1.split(".")]
        versions2 = [int(v) for v in version2.split(".")]
        for i in range(max(len(versions1), len(versions2))):
            v1 = versions1[i] if i < len(versions1) else 0
            v2 = versions2[i] if i < len(versions2) else 0
            if v1 > v2:
                return 1
            elif v1 < v2:
                return -1
        return 0


ob1 = Solution()

# Get package version
def geneutils_version():
    url = "https://pypi.org/project/geneutils/"
    source = requests.get(url)
    html_content = source.text
    soup = BeautifulSoup(html_content, "html.parser")
    company = soup.find("h1")
    vcheck = ob1.compareVersion(
        company.string.strip().split(" ")[-1],
        pkg_resources.get_distribution("geneutils").version,
    )
    if vcheck == 1:
        print(
            "\n"
            + "========================================================================="
        )
        print(
            "Current version of geneutils is {} upgrade to lastest version: {}".format(
                pkg_resources.get_distribution("geneutils").version,
                company.string.strip().split(" ")[-1],
            )
        )
        print(
            "========================================================================="
        )
    elif vcheck == -1:
        print(
            "\n"
            + "========================================================================="
        )
        print(
            "Possibly running staging code {} compared to pypi release {}".format(
                pkg_resources.get_distribution("geneutils").version,
                company.string.strip().split(" ")[-1],
            )
        )
        print(
            "========================================================================="
        )


geneutils_version()


accession_list = []
row_list = []
acession_annotations = []

# Save default NCBI email and API Key
def genauth():
    try:
        email_id = input("Enter your NCBI email address:  ")
        api_key = getpass.getpass("Enter your NCBI API Key:  ")
        if len(api_key) == 0:
            api_key = "pass"
    except Exception as e:
        print(e)

    data = {
        "api_key": api_key.strip(),
        "email": email_id.strip(),
    }

    with open(os.path.join(expanduser("~"), "geneauth.json"), "w") as outfile:
        json.dump(data, outfile)
    if len(data["api_key"]) > 0 and len(data["email"]) > 0:
        print("Credentials created with : email id & api key")


def genauth_from_parser(args):
    genauth()


# Merge lists together
def merge(lst1, lst2):
    return [a + b for (a, b) in zip(lst1, lst2)]


# Get annotate hit table
def accession_parse(path, db, email):
    """This script is intended for annotation of blast results saved in Hit Table CSV format.
       The output is an annotated CSV file "*_annotated.csv" with the following columns added:
       Record name, Species, Date of update, Reference, Full taxonomy

    [description]

    Arguments:
        path {[type]} -- [Pathway to csv-formatted Hit Table file with blast results]
        db {[type]} -- ["For the output of nucleotide blast or tblastn, use 'n'. "
                         "For the output of protein blast or blastx, use 'p'."]
        email {[type]} -- [NCBI email]
    """
    # drop empty rows
    df = pd.read_csv(path, header=None, usecols=[1])
    print("\n" + f"Total accession id(s): : {df.shape[0]}")
    df.dropna(inplace=True)
    df.drop_duplicates(keep="first", inplace=True)
    print(f"Total unique accession id(s): {df.shape[0]}")
    df.to_csv(path.split(".")[0] + "_deduped.csv")

    # get existing credentials
    if os.path.exists(os.path.join(expanduser("~"), "geneauth.json")):
        with open(os.path.join(expanduser("~"), "geneauth.json")) as json_file:
            cred = json.load(json_file)
        email_id = cred["email"].strip()
        api_key = cred["api_key"].strip()
    else:
        email_id = None
        api_key = None

    n = 50
    open(path.split(".")[0] + "_annotated.csv", "w")
    with open(path.split(".")[0] + "_deduped.csv") as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            accession_list.append(row[1])
            row_list.append(row)
    accession_chunk = [
        accession_list[i * n : (i + 1) * n]
        for i in range((len(accession_list) + n - 1) // n)
    ]
    i = 1
    length_accession_chunk = len(accession_chunk)
    if email is not None:
        Entrez.email = email
    elif email is None and email_id is not None:
        Entrez.email = email_id
    if api_key != "pass" and api_key is not None:
        Entrez.api_key = api_key
    if email is None and email_id is not None and api_key != "pass":
        print(f"Using both email and api_key: {email_id}")
    elif email is not None:
        print("Using email id {}".format(email))
    elif email is None and email_id is None:
        sys.exit("Either pass an email or try geneutils init")
    for c, subsets in enumerate(accession_chunk, 1):
        try:
            idlist = ",".join(subsets)
            # print(f'Processing {len(subsets)} accession id(s): ')
            dbdict = {"n": "nuccore", "p": "protein"}
            handle = Entrez.efetch(db=dbdict[db], id=idlist, retmode="xml")
            data = Entrez.parse(handle)
            for items in data:
                # break
                try:
                    version = items["GBSeq_accession-version"]
                    definition = items["GBSeq_definition"]
                    organism = items["GBSeq_organism"]
                    taxonomy = items["GBSeq_taxonomy"]
                    date = items["GBSeq_update-date"]
                    ref = ", ".join(items["GBSeq_references"][0]["GBReference_authors"])
                    with open(path.split(".")[0] + "_annotated.csv", "a") as csvfile:
                        writer = csv.writer(csvfile, delimiter=",", lineterminator="\n")
                        writer.writerow(
                            [version, definition, organism, taxonomy, ref, date]
                        )
                    csvfile.close()
                except Exception as e:
                    version = items["GBSeq_accession-version"]
                    definition = items["GBSeq_definition"]
                    organism = items["GBSeq_organism"]
                    taxonomy = items["GBSeq_taxonomy"]
                    date = items["GBSeq_update-date"]
                    ref = "No references found"
                    with open(path.split(".")[0] + "_annotated.csv", "a") as csvfile:
                        writer = csv.writer(csvfile, delimiter=",", lineterminator="\n")
                        writer.writerow(
                            [version, definition, organism, taxonomy, ref, date]
                        )
                    csvfile.close()
                except Exception:
                    print(f"Encountered exception with subset {c}")
                print(
                    f"Processed a total of {i} records: Subset {c} of {length_accession_chunk}",
                    end="\r",
                )
                i = i + 1
            handle.close()
        except Exception:
            pass
        except (KeyboardInterrupt, SystemExit) as e:
            sys.exit("Program escaped by User")

    print("")
    print("Annotated files written to {}".format(path.split(".")[0] + "_annotated.csv"))


def accession_parse_from_parser(args):
    accession_parse(path=args.path, db=args.db, email=args.email)


def main(args=None):
    parser = argparse.ArgumentParser(
        description="CLI and utilities for Genetic analysis and database interface"
    )
    subparsers = parser.add_subparsers()

    parser_genauth = subparsers.add_parser(
        "init", help="Setup NCBI email and API key as credentials"
    )
    parser_genauth.set_defaults(func=genauth_from_parser)

    parser_accession_parse = subparsers.add_parser(
        "blasthit", help="Annotation of blast results saved in Hit Table CSV format"
    )
    parser_accession_parse.add_argument(
        "--path", help="Pathway to csv-formatted Hit Table file with blast results"
    )
    parser_accession_parse.add_argument(
        "--db",
        choices=["n", "p"],
        help="For the output of nucleotide blast or tblastn, use 'n'. "
        "For the output of protein blast or blastx, use 'p'.",
    )
    optional_named = parser_accession_parse.add_argument_group(
        "Optional named arguments"
    )
    optional_named.add_argument(
        "--email", help="NCBI email else from credentials", default=None
    )
    parser_accession_parse.set_defaults(func=accession_parse_from_parser)

    args = parser.parse_args()

    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")
    func(args)


if __name__ == "__main__":
    main()
