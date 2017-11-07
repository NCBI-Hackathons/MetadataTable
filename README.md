# MetadataTable
Work on and around an interactive metadata table for sequence submission and modification

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