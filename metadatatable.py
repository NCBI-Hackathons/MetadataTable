# -*- coding: utf-8 -*-

import io
import sys
import csv
import time
import queue
import argparse
import threading
from Bio import Entrez
from metaparse import parse, get_experiment_xml_string
from metatestxml import TEST_XML

# ============================================================================ #
USAGE = """Usage examples:
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml -ot parsed.tsv -t biomol_transcript[properties] OR study_type_transcriptome_analysis[properties] OR strategy_rna_seq[properties] OR strategy_FL_cDNA[properties]
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml -ot parsed.tsv -t (staphylococcus aureus[Title]) AND blood[Text Word]
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml -t (staphylococcus aureus[Title]) AND blood[Text Word]
python metadatatable.py -e myemail@ncbi.gov -ot parsed.tsv -t strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov -t strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov -i myown.xml -ot parsed.tsv
python metadatatable.py -e myemail@ncbi.gov -i myown.xml -ot parsed.tsv -f
python metadatatable.py -e myemail@ncbi.gov -i myown.xml
"""

parser = argparse.ArgumentParser("MetadataTable", description=USAGE,
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-ox", "--output_xml",
                    type=str,
                    help="Path for saving xml file")

parser.add_argument("-ot", "--output_tsv",
                    type=str,
                    help="Path for saving parsed file")

parser.add_argument("-i", "--input_xml",
                    type=str,
                    help="Path to input file")

parser.add_argument("-e", "--email",
                    type=str,
                    help="Let NCBI know who you are (required if using -t)")

parser.add_argument("-t", "--term",
                    nargs="+",
                    type=str,
                    help="Query terms")

parser.add_argument("-u", "--unlimited",
                    action="store_true",
                    help="Retrieve unlimited records")

parser.add_argument("-c", "--case",
                    type=str,
                    default="rnaseq",
                    choices=("rnaseq", "source"),
                    help="Select which builtin case to use")

parser.add_argument("-f", "--full",
                    action="store_true",
                    help="Whether to output full table (Only for builtin template)")

parser.add_argument("-x", "--xpath",
                    type=str,
                    help="Path to a csv file for xpath query")

# ============================================================================ #
BUILTIN_XPATH = {
    "rnaseq": [
        ("SRA run accession", ("//RUN_SET/RUN", "accession"), True),
        ("SRA experiment accession", ("//EXPERIMENT", "accession", "accession"), True),
        ("Biosample accession (1-to-1 with SRA sample accession when both exist)", ("//IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']", None), True),
        ("Tissue", ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue']/VALUE|//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue source']/VALUE|//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='OrganismPart']/VALUE", None), True),
        ("Strain", ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='strain']/VALUE", None), True),
        ("Development Stage", ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='developmental stage']/VALUE|//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='DevelopmentalStage']/VALUE", None), True),
        ("SRA project accession", ("//STUDY_REF", "accession"), True),
        ("Base count of run", ("//RUN", "total_bases"), True),
        ("Paired-end flag", ("//PAIRED","*"), True),
        ("Spot count of run", ("//RUN", "total_spots"), True),
        ("Platform (eg Illumina)", ("//PLATFORM/*/*|//PLATFORM/*", None), True),
        ("SRA sample accession", ("//SAMPLE","accession"), False),
        # ("Taxid", ("todo",), False),
        # ("Library source", ("todo",), False),
        # ("Cell line", ("todo",), False),
        # ("Sample title", ("todo",), False),
        # ("Source Provider", ("todo",), False),
        # ("Study description", ("todo",), False),
    ],
    "source": [
    ],
}

# ============================================================================ #
ESEARCH_BATCH = 100000  # max=100,000
EFETCH_BATCH = 1000  # max=10,000
ESEARCH_MAX = ESEARCH_BATCH


def Wait(lastTime):  # Please do not post more than three URL requests per second.
    while time.time() - lastTime < 0.33:
        time.sleep(0.11)
    return time.time()


def GetIdList(**kwargs):
    idList = []
    # Get total count
    handle = Entrez.esearch(db="sra", retmax=ESEARCH_BATCH, **kwargs)
    result = Entrez.read(handle)
    handle.close()
    idList.extend(result["IdList"])
    total = int(result["Count"])
    # Get all idList
    retstart = ESEARCH_BATCH
    lastTime = time.time()
    while retstart < total:
        handle = Entrez.esearch(db="sra", retmax=ESEARCH_BATCH, retstart=retstart, **kwargs)
        result = Entrez.read(handle)
        handle.close()
        idList.extend(result["IdList"])
        retstart += ESEARCH_BATCH
        lastTime = Wait(lastTime)
    return idList


def GetRecords(idList):
    total = len(idList)
    retstart = 0
    lastTime = time.time()
    while retstart < total:
        handle = Entrez.efetch(db="sra", id=idList, retmax=EFETCH_BATCH, retstart=retstart)
        data = io.TextIOWrapper(handle.detach(), encoding="utf-8")
        yield data.read()
        retstart += EFETCH_BATCH
        lastTime = Wait(lastTime)


# ============================================================================ #
class TaskDone(object):
    pass


class BlackHole(object):
    def write(self, value):
        pass


def Process(q, output_xml, output_tsv, names, queries):
    fx = ft = BlackHole()
    if output_xml:
        fx = open(output_xml, "w", encoding="utf-8")
        fx.write('<?xml version="1.0" encoding="UTF-8"?>\n<EXPERIMENT_PACKAGE_SET>\n')
    if output_tsv:
        ft = open(output_tsv, "w", encoding="utf-8")
        ft.write("\t".join(names))
        ft.write("\n")
    if output_xml is None and output_tsv is None:
        ft = sys.stdout
    while True:
        item = q.get()
        if item is TaskDone:
            q.task_done()
            break
        fx.write("".join(get_experiment_xml_string(item)))
        for row in parse(item, queries):
            ft.write("\t".join(row))
            ft.write("\n")
            ft.flush()
        q.task_done()
    if output_xml:
        fx.write("</EXPERIMENT_PACKAGE_SET>\n")
        fx.close()
    if output_tsv:
        ft.close()


# ============================================================================ #
if __name__ == "__main__":
    P, _ = parser.parse_known_args()

    # print("P.output_xml", P.output_xml)
    # print("P.output_tsv", P.output_tsv)
    # print("P.input_xml", P.input_xml)
    # print("P.email", P.email)
    # print("P.term", P.term)
    # print("P.unlimited", P.unlimited)
    # print("P.case", P.case)
    # print("P.full", P.full)
    # print("P.xpath", P.xpath)
    # exit()

    # Only one input mode is allowed to avoid confusion
    if P.input_xml and P.term:
        sys.stderr.write("Please use either a xml file or search terms as input")
        exit(1)
    if P.term and not P.email:
        sys.stderr.write("Email is required for querying Entrez")
        exit(1)

    # Set the field names and queries
    if P.xpath:
        NAMES = []
        QUERIES = []
        with open(P.xpath, newline="", encoding="utf-8") as csvfile:
            for row in csv.reader(csvfile):
                NAMES.append(row[0])
                QUERIES.append((row[1], row[2] if len(row) == 3 and row[2] else None))
        try:
            for _ in parse(TEST_XML, QUERIES):
                pass
        except Exception as e:
            sys.stderr.write("Xpath test failed. Please check the xpath file.")
            exit(1)
    else:
        fields = BUILTIN_XPATH[P.case]
        if P.full:
            NAMES = [row[0] for row in fields if row[2]]
            QUERIES = [row[1] for row in fields if row[2]]
        else:
            NAMES = [row[0] for row in fields]
            QUERIES = [row[1] for row in fields]

    # Input Mode 1. Read and process xml file, save parsed data to a file or print out parsed data
    if P.input_xml:
        if P.output_tsv:
            fo = open(P.output_tsv, "w", encoding="utf-8")
            fo.write("\t".join(NAMES))
            fo.write("\n")
        else:
            fo = sys.stdout
        with open(P.input_xml, "r", encoding="utf-8") as fi:
            for result in parse(fi, QUERIES):
                fo.write("\t".join(result))
                fo.write("\n")
                fo.flush()
        if P.output_tsv:
            fo.close()

    # Input Mode 2. Query Entrez and parse on the fly
    elif P.term:
        Entrez.email = P.email
        Entrez.tool = "MetadataTable"

        # Retrieve all ids
        ids = GetIdList(term=" ".join(P.term))
        if len(ids) > ESEARCH_MAX:
            sys.stderr.write("Query returned too many results (%s). Please consider refine you search or use -u option" % len(ids))
            exit(1)

        # Process the downloaded records in a separate thread
        q = queue.Queue()
        t = threading.Thread(target=Process, args=(q, P.output_xml, P.output_tsv, NAMES, QUERIES), daemon=True)
        t.start()

        # Retrieve all records and put them in the queue
        for d in GetRecords(ids):
            q.put(d)
        q.put(TaskDone)

        # Wait till all records are processed
        q.join()

    exit(0)
