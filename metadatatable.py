# -*- coding: utf-8 -*-

import io
import sys
import time
import queue
import argparse
import threading
from Bio import Entrez

from parse import parse, generate_xml_string

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

parser.add_argument("-e", "--email",
                    type=str,
                    help="Let NCBI know who you are (required)")

parser.add_argument("-ox", "--output_xml",
                    type=str,
                    help="Path for saving xml file")

parser.add_argument("-ot", "--output_tsv",
                    type=str,
                    help="Path for saving parsed file")

parser.add_argument("-i", "--input_xml",
                    type=str,
                    help="Path to input file")

parser.add_argument("-f", "--full",
                    action="store_true",
                    help="Whether to output full table")

parser.add_argument("-t", "--term",
                    nargs="+",
                    type=str,
                    help="Query terms")

# ============================================================================ #
ESEARCH_BATCH = 100000  # max=100,000
EFETCH_BATCH = 1000  # max=10,000


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
        fx.write(generate_xml_string(item))
        for row in parse(item, queries):
            ft.write("\t".join(row))
            ft.write("\n")
        q.task_done()
    if output_xml:
        fx.write("</EXPERIMENT_PACKAGE_SET>\n")
        fx.close()
    if output_tsv:
        ft.close()


# ============================================================================ #
if __name__ == "__main__":
    # parsed, _ = parser.parse_known_args("-e yourmail@here.com -ox sra1000.xml txid112509[Organism:exp]".split())
    P, _ = parser.parse_known_args()

    # print(P.email)
    # print(P.input_xml)
    # print(P.output_xml)
    # print(P.output_tsv)
    # print(P.term)
    # print(P.full)
    # exit()

    if P.input_xml and P.term:
        sys.stderr.write("Please use either a xml file or search terms as input")
        exit()

    FIELDS = [  # (Column name, query, default)
        ("SRA run accession", (".//RUN_SET/RUN", "accession"), True),
        ("SRA experiment accession", ("../../EXPERIMENT", "accession"), True),
        ("Biosample accession (1-to-1 with SRA sample accession when both exist)", "../../IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']", True),
        #("Tissue", ("todo",), True),
        #("Strain", ("todo",), True),
        #("Developmental stage", ("todo",), True),
        #("SRA project accession", ("todo",), True),
        #("Base count of run", ("todo",), True),
        #("Spot count of run", ("todo",), True),
        #("Paired-end flag", ("todo",), True),
        #("Platform (eg Illumina)", ("todo",), True),
        #("SRA sample accession", ("todo",), False),
        #("Taxid", ("todo",), False),
        #("Library source", ("todo",), False),
        #("Cell line", ("todo",), False),
        #("Sample title", ("todo",), False),
        #("Source Provider", ("todo",), False),
        #("Study description", ("todo",), False),
    ]
    if not P.full:
        NAMES = [row[0] for row in FIELDS]
        QUERIES = [row[1] for row in FIELDS]
    else:
        NAMES = [row[0] for row in FIELDS if row[2]]
        QUERIES = [row[1] for row in FIELDS if row[2]]

    # Input Mode 1. Read and process xml file
    if P.input_xml:
        # Save parsed data to a file or just print out parsed data
        if P.output_tsv:
            fo = open(P.output_tsv, "w", encoding="utf-8")
            fo.write("\t".join(NAMES))
            fo.write("\n")
        else:
            fo = sys.stdout
        for result in parse(open(P.input_xml, "r", encoding="utf-8"), QUERIES):
            fo.write("\t".join(result))
            fo.write("\n")
        if P.output_tsv:
            fo.close()
    # Input Mode 2. Query Entrez and parse on the fly
    elif P.term:
        Entrez.email = P.email
        Entrez.tool = "MetadataTable"
        q = queue.Queue()
        # Process the downloaded records in a separate thread
        t = threading.Thread(target=Process, args=(q, P.output_xml, P.output_tsv, NAMES, QUERIES))
        t.start()
        # Retrieve all ids
        ids = GetIdList(term=" ".join(P.term))[:10]  # TODO
        # Retrieve all records and put them in the queue
        for d in GetRecords(ids):
            q.put(d)
        q.put(TaskDone)
        # Wait till all records are processed
        q.join()