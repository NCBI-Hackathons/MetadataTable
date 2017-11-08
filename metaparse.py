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
    corrected_experiments = []
    outer_xml = re.split("</?EXPERIMENT_PACKAGE_SET>", xml_content)
    experiments = re.split("<EXPERIMENT_PACKAGE>", outer_xml[1].strip())
    for experiment in experiments[1:]:
        corrected_experiments.append("<EXPERIMENT_PACKAGE>" + experiment)
    return corrected_experiments


def parse(file_or_string_to_parse, xpath_list):
    """returns a list of values requested from
       the xpath_list for each run in the xml string or file"""
    xml_content = _return_file_content(file_or_string_to_parse)
    experiments = get_experiment_xml_string(xml_content)
    for experiment in experiments:
        tree = etree.fromstring(experiment) #loads each experiment into a tree
        runs = tree.findall(".//RUN_SET/RUN") #obtains all runs from the experiment
        list_of_values = []
        for i, run in enumerate(runs):
            list_of_values = []
            print(i)
            for field in xpath_list:
                run_relative_xpath = "./../.." + field[0]
                el = run.find(run_relative_xpath)
                if el is None:
                    list_of_values.append("NA")
                elif field[1]:
                        list_of_values.append(el.attrib.get(field[1], "NA"))
                else:
                        list_of_values.append(el.text)

            yield list_of_values #yields a list for each run
