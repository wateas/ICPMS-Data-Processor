# ICPMS-Data-Processor
This project utilizes python to combine and process dataframes sourced from excel spreadsheets of ICP-MS data that have been exported from Agilent's Masshunter software.  A series of excel spreadsheets demonstrating the transformations, a customized reporting/quantification limit chart, and a table of extracted method blanks are returned.

## Getting Started

Several excel workbooks have been included that can be used as test data.  

The csv file "quant_limit_chart_ICPMS" should be in the same directory as the source data.  The user should select appropriate values that correspond to their instrument/application.

#### Prerequisites

A python environment file (.yml) has been included.

Use to following in anaconda prompt to recreate the environment:

`% conda env create -n conda-env-name -f /path/to/environment.yml`

## Future Versions

Future versions will further automate QC validation and data processing.
