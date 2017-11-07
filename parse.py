from lxml import etree
import re

def _return_file_content(file_or_string_to_parse):
    if isinstance(file_or_string_to_parse, str):
        return file_or_string_to_parse
    else:
        return file_or_string_to_parse.read()


def parse(file_or_string_to_parse):
    data_array = []
    exp_dict = {}
    xml_content = _return_file_content(file_or_string_to_parse)
    outer_xml = re.split("</?EXPERIMENT_PACKAGE_SET>", xml_content)
    experiments = re.split("<EXPERIMENT_PACKAGE>", outer_xml[1].strip())
    for experiment in experiments[1:]:
        experiment = "<EXPERIMENT_PACKAGE>" + experiment
        tree = etree.fromstring(experiment)
        runs = tree.findall(".//RUN_SET/RUN")
        list_of_values = []
        for run in runs:
            for field in [("./../..//STUDY_REF", "accession"),
                          (".","total_spots")]:
                if isinstance(field, tuple):
                    list_of_values.append(run.find(field[0]).attrib.get(field[1], "NA"))
                else:
                    list_of_values.append(run.find(field).text)
            yield list_of_values

sra = open("sra1000.xml", "r")
parse_func = (parse(sra))
print(next(parse_func))
