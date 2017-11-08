from lxml import etree
import re

def _return_file_content(file_or_string_to_parse):
    if isinstance(file_or_string_to_parse, str):
        return file_or_string_to_parse
    else:
        return file_or_string_to_parse.read()

def get_experiment_xml_string(xml_content):
    outer_xml = re.split("</?EXPERIMENT_PACKAGE_SET>", xml_content)
    experiments = re.split("<EXPERIMENT_PACKAGE>", outer_xml[1].strip())
    return experiments[1:]


def parse(file_or_string_to_parse, xpath_list):
    data_array = []
    exp_dict = {}
    xml_content = _return_file_content(file_or_string_to_parse)
    experiments = get_experiment_xml_string(xml_content)
    for experiment in experiments:
        experiment = "<EXPERIMENT_PACKAGE>" + experiment
        tree = etree.fromstring(experiment)
        runs = tree.findall(".//RUN_SET/RUN")
        list_of_values = []
        for run in runs:
            for field in xpath_list:
                if isinstance(field, tuple):
                    list_of_values.append(run.find(field[0]).attrib.get(field[1], "NA"))
                else:
                    list_of_values.append(run.find(field).text)
            yield list_of_values

sra = open("sra1000.xml", "r")
parse_func = (parse(sra))
print(next(parse_func))
