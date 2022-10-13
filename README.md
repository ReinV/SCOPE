# SCOPE
<!-- &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.3 Workflow](#13Workflow)   -->
[1. Introduction](#1-Introduction)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.1 What is SCOPE?](#11-What-is-SCOPE)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.2 How does it work?](#12-How-does-it-work?)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.3 What can it be used for?](#13-What-can-it-be-used-for)  
[2. Prepare SCOPE](#2-Prepare-SCOPE)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.1 Pull the git repo](#21-Pull-the-git-repo)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.2 Download the required files](#22-Download-the-required-files)  
[3. Usage](#3-Usage)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1 Create query file](#31-Create-query-file)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2 Search for all publications using the query file](#32-Search-for-all-publications-using-the-query-file)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.3 Summarize the query results: make a table](#33-Summarize-the-query-results:-make-a-table)  
[4. Further reading](#4-Further-reading)  

## Introduction 

## 1.1 What is SCOPE?
SCOPE is a literature review tool for detecting chemical patterns. Chemicals resulting from a query search can be visualized in an interactive density plot (mass over logP). The more a chemical is mentioned, the brighter the spot. More information is shown in a table by hovering over a spot. The user can easily switch between plots, and add blur and saturation. This allows for comparison of mass/logP distributions between different query searches. 

## 1.2 How does it work?
Chemicals in publications on the [Europe PMC site](https://europepmc.org/) are tagged with a CheBI identifier.  Several chemical properties can be retrieved from the [ChEBI Ontology](https://www.ebi.ac.uk/chebi/) using a ChEBI identifier as a key. Log *P* and log *S* values are predicted using the AlogPS3.0 model from the OCHEM website [ochem.eu site](https://ochem.eu) with SMILES (notification for chemical structures, also retrieved from the ontology). SCOPE connects ChEBI identifiers resulting from a query results with chemical properties of interest.

<!-- -- workflow image --  -->

## 1.3 What can it be used for?
- Compare chemical detection methods
- Compare vendors
- Compare medicines/diseases

<!-- 
This toolbox contains python scripts that follow the [EuropePMC2ChEBI KNIME workflow](https://github.com/magnuspalmblad/EuropePMC2ChEBI). In this workflow, annotated chemicals are collected through literature searches via the [Europe PMC site](https://europepmc.org/). Properties of these chemicals can be visualized in interactive plots using the [Bokeh library](https://bokeh.pydata.org). This allows comparison between the data of different query searches, for example detecting bias in different analytical chemical techniques.

# Plotting chemicals
For plotting chemicals we use certain properties such as mass and logP. These do not come from the search itself but have been retrieved from several sources. We use the ChEBI Ontology to retrieve all ChEBI chemicals, their mass values, and SMILES. Additionally, log *P* and log *S* values are predicted using the AlogPS3.0 model from the OCHEM website [ochem.eu site](https://ochem.eu). -->

# 2. Prepare SCOPE

## 2.1 Pull the git repo

using ssh 

<pre><code>git pull git@github.com:ReinV/SCOPE.git</code></pre>

or https

<pre><code>git pull https://github.com/ReinV/SCOPE.git</code></pre>

## 2.2 Download the required files from our [OSF project](https://osf.io/pvwu2/)

<pre><code>python download_files.py</code></pre>

In these files, ChEBI identifiers are linked to properties of interest (for plotting).

Everything should now be in place for using SCOPE.

# 3. Usage

Execute "python \<script name> -h" to get a description of the required arguments.

## 3.1 Create query file

Put the search queries (at least one) in a text file in the following structure: 

<pre><code>&lt;output-tag-1&gt;, &lt;query-1&gt;
&lt;output-tag-2&gt;, &lt;query-2&gt;
&lt;output-tag-3&gt;, &lt;query-3&gt;</code></pre>

The search query should be written with the syntax used on the Europe PMC site. For example:
> METHODS:"Nuclear Magnetic Resonance" OR METHODS:NMR OR METHODS:"NMR spectrometry" OR METHODS:"nuclear magnetic resonance spectrometry" OR METHODS:"NMR spectroscopy" OR METHODS:"nuclear magnetic resonance (NMR) spectroscopy"

We provide one example input file in the "queries" folder. 

## 3.2 Search for all publications using the query file

<pre><code>python search_query.py -i &lt;path-to-input-file&gt;</code></pre>

The results will be stored in the "results" folder with the output tag as output name. 

Warning: this may take up to several hours if there are many search hits. It is recommended to start with a more specific search when trying out SCOPE for the first time.

## 3.3 Summarize the query results: make a table

Use the results folder as input

<pre><code>python make_table.py -i results -t folder</code></pre>

or a specific file

<pre><code>python make_table.py -i &lt;path-to-result-file&gt; -t file</code></pre>

## 3.4 Plot the query results

<pre><code>python visualize_multiplot.py -i tables -o &lt;plot-name&gt;</code></pre>

This plot will be saved in the "plots" folder.

# 4 Further reading


