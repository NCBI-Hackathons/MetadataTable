# -*- coding: utf-8 -*-

import io
import re
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

# parser.add_argument("-c", "--case",
#                     type=str,
#                     default="rnaseq",
#                     choices=("rnaseq", "source"),
#                     help="Select which builtin case to use")

parser.add_argument("-f", "--full",
                    action="store_true",
                    help="Whether to output full table (Only for builtin template)")

parser.add_argument("-x", "--xpath",
                    type=str,
                    help="Path to a csv file for xpath query")

# ============================================================================ #
SHOW = True
HIDE = False
TEXT = False

BUILTIN_XPATH = {
    "rnaseq": [
        # Accessions
        (SHOW, "Project Accession", (("//STUDY_REF/IDENTIFIERS/PRIMARY_ID", TEXT),)),  # ("//STUDY_REF", "accession")
        (HIDE, "Project Accession (Secondary)", (
            ("//STUDY_REF/IDENTIFIERS/EXTERNAL_ID[@namespace='BioProject']", TEXT),
            ("//STUDY/IDENTIFIERS/EXTERNAL_ID[@namespace='BioProject']", TEXT),
        )),
        (SHOW, "Sample Accession", (("//SAMPLE/IDENTIFIERS/PRIMARY_ID", TEXT),)),  # ("//SAMPLE", "accession")
        (HIDE, "Sample Accession (Secondary)", (("//SAMPLE/IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']", TEXT),)),
        (SHOW, "Experiment Accession", (("//EXPERIMENT/IDENTIFIERS/PRIMARY_ID", TEXT),)),  # ("//EXPERIMENT", "accession")
        (SHOW, "Run Accession", (("//RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID", TEXT),)),  # ("//RUN_SET/RUN", "accession")

        # Design
        (SHOW, "Platform", (
            ("//PLATFORM/*/*", TEXT),
            ("//PLATFORM/*", TEXT),
        )),
        (SHOW, "Layout", (
            ("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_LAYOUT/PAIRED", ("Paired", "")),
            ("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_LAYOUT/SINGLE", ("Single", "")),
        )),
        (HIDE, "Library Name", (("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_NAME", TEXT),)),
        (HIDE, "Library Strategy", (("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_STRATEGY", TEXT),)),
        (HIDE, "Library Source", (("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_SOURCE", TEXT),)),
        (HIDE, "Library Selection", (("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_SELECTION", TEXT),)),
        (HIDE, "Library Construction", (("//DESIGN//LIBRARY_DESCRIPTOR//LIBRARY_CONSTRUCTION_PROTOCOL", TEXT),)),

        # Study
        (HIDE, "Study Title", (("//STUDY//DESCRIPTOR//STUDY_TITLE", TEXT),)),
        (HIDE, "Study Type", (("//STUDY//DESCRIPTOR//STUDY_TYPE", "existing_study_type"),)),
        (SHOW, "Study Abstract", (("//STUDY//DESCRIPTOR//STUDY_ABSTRACT", TEXT),)),

        # Sample
        (HIDE, "Sample Title", (
            ("//SAMPLE//TITLE", TEXT),
            ("//RUN//Pool/Member", "sample_title"),
        )),
        (HIDE, "Taxon Id", (("//SAMPLE//SAMPLE_NAME//TAXON_ID", TEXT),)),
        (HIDE, "Scientific Name", (("//SAMPLE//SAMPLE_NAME//SCIENTIFIC_NAME", TEXT),)),
        (SHOW, "Strain/Cultivar", (
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='strain']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='genetic background']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='cultivar']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='ecotype']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='varieties']/VALUE", TEXT),
        )),
        (HIDE, "Genotype", (
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='genotype']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='genotype/variation']/VALUE", TEXT),
        )),
        (HIDE, "Phenotype", (("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='phenotype']/VALUE", TEXT),)),
        (SHOW, "Tissue", (
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='tissue source']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='OrganismPart']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='organism part']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='source_name']/VALUE", TEXT),
        )),
        (SHOW, "Development Stage", (
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='developmental stage']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='development stage']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='DevelopmentalStage']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='dev_stage']/VALUE", TEXT),
            ("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='Stage']/VALUE", TEXT),
        )),
        (HIDE, "Age Value", (("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='age']/VALUE", TEXT),)),
        (HIDE, "Age Unit", (("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='age']/UNITS", TEXT),)),
        (HIDE, "Treatment", (("//SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='treatment']/UNITS", TEXT),)),

        # Run
        (SHOW, "Read Length", (("//RUN/Statistics/Read", "average"),)),
        (SHOW, "Total Bases", (("//RUN", "total_bases"),)),
        (SHOW, "Total Spots", (("//RUN", "total_spots"),)),

        # Misc
        (HIDE, "ArrayExpress ID", (("//STUDY/STUDY_ATTRIBUTES/STUDY_ATTRIBUTE[TAG='ArrayExpress']/VALUE", TEXT),)),
        (SHOW, "Published", (("//RUN_SET/RUN", "published"),)),

        # "Cell line",
        # "Source Provider",
    ],
    "source": [],  # NotImplemented
}

RE_RUN = re.compile(r'(RUN\b)|(RUN\/)')


def TransformXpath(xpath):
    run_container = RE_RUN.search(xpath)
    if not run_container:
        return "./../.." + xpath
    else:
        return "." + xpath[run_container.end():]


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

    def flush(self):
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
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

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
    # sys.exit(1)

    # Only one input mode is allowed to avoid confusion
    if P.input_xml and P.term:
        sys.stderr.write("Please use either a xml file or search terms as input")
        sys.exit(1)
    # Email is required for e-utils
    if P.term and not P.email:
        sys.stderr.write("Email is required for querying Entrez")
        sys.exit(1)

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
            sys.exit(1)
    else:
        # fields = BUILTIN_XPATH[P.case]
        fields = BUILTIN_XPATH["rnaseq"]
        if P.full:
            NAMES = [row[1] for row in fields]
            QUERIES = [row[2] for row in fields]
        else:
            NAMES = [row[1] for row in fields if row[0]]
            QUERIES = [row[2] for row in fields if row[0]]
        QUERIES = [[(TransformXpath(q[0]), q[1]) for q in query] for query in QUERIES]

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
        if len(ids) > ESEARCH_MAX and not P.unlimited:
            sys.stderr.write("Query returned too many results (%s). Please consider refine you search or use -u option" % len(ids))
            sys.exit(1)

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

    sys.exit(0)
