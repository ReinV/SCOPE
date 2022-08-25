# Toolbox
This toolbox contains python scripts that follow the [EuropePMC2ChEBI KNIME workflow](https://github.com/magnuspalmblad/EuropePMC2ChEBI). In this workflow, annotated chemicals are collected through literature searches via the [Europe PMC site](https://europepmc.org/). Properties of these chemicals can be visualized in interactive plots using the [Bokeh library](https://bokeh.pydata.org). This allows comparison between the data of different query searches, for example detecting bias in different analytical chemical techniques.

# Plotting chemicals
For plotting chemicals we use certain properties such as mass and logP. These do not come from the search itself but have been retrieved from several sources. We use the ChEBI Ontology to retrieve all ChEBI chemicals, their mass values, and SMILES. Additionally, log *P* and log *S* values are predicted using the AlogPS3.0 model from the OCHEM website [ochem.eu site](https://ochem.eu).

# Publication 
In progress ...

<!-- # Workflow -->
<!-- ![Workflow scheme](workflow_scheme_new_adjusted.png) -->

# Prepare SCOPE

1. Pull the git repo

using ssh 

<pre><code>git pull git@github.com:ReinV/SCOPE.git</code></pre>

or https

<pre><code>git pull https://github.com/ReinV/SCOPE.git</code></pre>


2. Download the required files from our [OSF project](https://osf.io/pvwu2/). 

<pre><code>python download_files.py</code></pre>

In these files, ChEBI identifiers are linked to properties of interest (for plotting).

Everything should now be in place for using SCOPE.

# Usage

Execute "python \<script name> -h" to get a description of the required arguments.

1. Create query file

Put the search queries (at least one) in a text file in the following structure: 

<pre><code>&lt;output-tag-1&gt;, &lt;query-1&gt;
&lt;output-tag-2&gt;, &lt;query-2&gt;
&lt;output-tag-3&gt;, &lt;query-3&gt;</code></pre>

The search query should be written with the syntax used on the Europe PMC site. For example:
> METHODS:"Nuclear Magnetic Resonance" OR METHODS:NMR OR METHODS:"NMR spectrometry" OR METHODS:"nuclear magnetic resonance spectrometry" OR METHODS:"NMR spectroscopy" OR METHODS:"nuclear magnetic resonance (NMR) spectroscopy"

We provide one example input file in the "queries" folder. 

2. Search for all publications using the query file. 

<pre><code>python search_query.py -i &lt;path-to-input-file&gt;</code></pre>

The results will be stored in the "results" folder with the output tag as output name. 

Warning: this may take up to several hours if there are many search hits. It is recommended to start with a more specific search when trying out SCOPE for the first time.

3. Make a table. 

Use the results folder as input

<pre><code>python make_table.py -i results -t folder</code></pre>

or a specific file

<pre><code>python make_table.py -i &lt;path-to-result-file&gt; -t file</code></pre>

4. Plot the query results. 

<pre><code>python visualize_multiplot.py -i tables -o &lt;plot-name&gt;</code></pre>

This plot will be saved in the "plots" folder.


