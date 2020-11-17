# ICPMS-Data-Processor
Automates elemental data analysis.

The purpose of this script is to combine and process dataframes sourced from excel 
spreadsheets of ICP-MS data that have been exported from Agilent's Masshunter 
software.

The first several steps of processing ICP-MS or other elemental data is performed
and a series of excel spreadsheets demonstrating transformations, a 
reporting/quantification limit chart, and extracted method blanks are returned.

The csv file "quant_limit_chart_ICPMS" should be in the same directory as the source 
data.  The user should select appropriate values that correspond to their 
instrument/application.

Future versions will further automate QC validation and data processing.
