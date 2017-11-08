# MetadataTable
Command line resource to retrieve a subset of SRA metadata in table form.  Intended to be easily modifiable and flexible so that users can identify SRA projects of interest.

# What's the problem?
Researchers who seek to identify SRA projects of interest using more detailed criteria than those supported by Entrez often implement their own scripted or database solution to examine metadata of many studies.  This script allows researchers uninterested in scripting to retrieve a summary of SRA metadata for many studies.  For example, one seeks to browse all RNA-seq studies from a particular organism to identify those that include a particular tissue or developmental stage.  This particular example is preconfigured in the script.  

One technical obstacle to easily querying SRA is the heterogeneity of the data.  There are for example several ways to specify that a study contains RNA-Seq rather than WGS or exome data, and further several ways to identify that an SRA sample consists of a tissue of interest.  For some researchers, several terms are necesary to write a biologically useful query to retrieve ALL SRA entities of interest, rather than just those which satisfy a single specification.  For RNA-Seq studies and samples originating from particular tissues, the synonymous terms are preset.  Interested users can override the defaults or modify the relatively straightforward code in this repository.  

For users interested to identify terms of interest-- browsing the SAMPLE_ATTRIBUTE fields of selected SRA records is often a useful way to identify how studies are described.


# Usage
```
>python metadatatable.py --help
usage: MetadataTable [-h] -e EMAIL [-ox OUTPUT_XML] [-ot OUTPUT_TSV]
                     [-i INPUT_XML]
                     term [term ...]

Usage examples:
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml -ot parsed.tsv biomol_transcript[properties] OR study_type_transcriptome_analysis[properties] OR strategy_rna_seq[properties] OR strategy_FL_cDNA[properties]
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml -ot parsed.tsv strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov -ox raw.xml strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov -ot parsed.tsv strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov strategy_rna_seq[properties]
python metadatatable.py -e myemail@ncbi.gov -i myown.xml -ot parsed.tsv
python metadatatable.py -e myemail@ncbi.gov -i myown.xml

positional arguments:
  term                  Query terms

optional arguments:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL
                        Let NCBI know who you are
  -ox OUTPUT_XML, --output_xml OUTPUT_XML
                        Path for saving xml file
  -ot OUTPUT_TSV, --output_tsv OUTPUT_TSV
                        Path for saving parsed file
  -i INPUT_XML, --input_xml INPUT_XML
                        Path to input file
```
