# Toolbox
This toolbox contains python scripts that follow the [EuropePMC2ChEBI KNIME workflow](https://github.com/magnuspalmblad/EuropePMC2ChEBI). In this workflow, annotated chemicals are collected through literature searches via the [Europe PMC site](https://europepmc.org/). Properties of these chemicals can be visualized in interactive plots using the [Bokeh library](https://bokeh.pydata.org). This allows comparison between the data of different query searches, for example detecting bias in different analytical chemical techniques.

# Workflow
![Workflow scheme](workflow_scheme_new_adjusted.png)

# Preparing SCOPE

1. Run download_files.py. This script will download the required files from the [OSF project](https://osf.io/pvwu2/) and create subfolders under the current folder.

All files should now be in place for running SCOPE.

# Usage

Execute "python \<script name> -h" to get a description of the required arguments.

1. (Optional) Execute "python update_chebis.py". This script will check if the property files in the "files" folder can be updated with the latest ChEBI ontology.

2. Put the search queries (at least one) in a text file in the following structure: </br><output_tag_1>, <search query 1>  
<output_tag_2>, <search query 2>  
<output_tag_3>, <search query 3>  

The search query should be written with the syntax used on the Europe PMC site (see ... ). For example:
> METHODS:"Nuclear Magnetic Resonance" OR METHODS:NMR OR METHODS:"NMR spectrometry" OR METHODS:"nuclear magnetic resonance spectrometry" OR METHODS:"NMR spectroscopy" OR METHODS:"nuclear magnetic resonance (NMR) spectroscopy"

We provide one example text file in the "queries" folder.

4. Execute "python search_query.py -i queries/\<input text file>" to search for all publications. The results will be stored in the "results" folder with the output tag as output name. Warning: this may take up many hours if there are a lot of search hits!

5. Execute "python make_table.py -i results -t folder" to create tables in the "tables" folder for all results.

6. Execute "python visualize_multiplot.py -i tables -o <output name>" to create a plot using the tables in the "tables" folder. This plot will be saved in the "plots" folder.

# Using external sources to get plot properties
In the 'files' folder, ChEBI identifiers are linked to a certain property. These do not come from the search itself but can be looked up in the ChEBI Ontology. Additionally, log *P* and log *S* values are predicted using the AlogPS3.0 model from the OCHEM website [ochem.eu site](https://ochem.eu).
