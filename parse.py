from lxml import etree
import re

def _return_file_content(file_or_string_to_parse):
    """checks if input is a file object and returns a string if it is"""
    if isinstance(file_or_string_to_parse, str):
        return file_or_string_to_parse
    else:
        return file_or_string_to_parse.read()

def get_experiment_xml_string(xml_content):
    """"returns an xml string of experiments
        outside the experiment package set root"""
    outer_xml = re.split("</?EXPERIMENT_PACKAGE_SET>", xml_content)
    experiments = re.split("<EXPERIMENT_PACKAGE>", outer_xml[1].strip())
    return experiments[1:]


def parse(file_or_string_to_parse, xpath_list):
    """returns a list of values requested from
       the xpath_list for each run in the xml string or file"""
    data_array = []
    exp_dict = {}
    xml_content = _return_file_content(file_or_string_to_parse)
    experiments = get_experiment_xml_string(xml_content)
    for experiment in experiments:
        experiment = "<EXPERIMENT_PACKAGE>" + experiment
        tree = etree.fromstring(experiment) #loads each experiment into a tree
        runs = tree.findall(".//RUN_SET/RUN") #obtains all runs from the experiment
        list_of_values = []
        for run in runs:
            for field in xpath_list:
                if isinstance(field, tuple): #checks if the xpath is a tuple with an attribute
                    list_of_values.append(run.find(field[0]).attrib.get(field[1], "NA")) #obtains the attribute
                else:
                    list_of_values.append(run.find(field).text) #adds all values to a list
            yield list_of_values #yields a list for each run
