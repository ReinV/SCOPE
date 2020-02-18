#!/usr/bin/python
import time

import numpy as np
import pandas as pd
import argparse
from math import exp
from math import sqrt
from datetime import datetime
from os import listdir
import sys

from bokeh import events
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.models import CustomJS
from bokeh.models import HoverTool
from bokeh.models import CustomJS, ColumnDataSource, Slider, CheckboxGroup, RadioGroup, Button, MultiSelect
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import axial_to_cartesian
from bokeh.util.hex import cartesian_to_axial
from bokeh.layouts import column, row
from bokeh.palettes import Viridis256, Greys256

def import_table(file):
    '''
    This function imports the table as csv file and returns a dictionary with the ChEBI ID's as keys and a dictionary with their poperties as values.
    '''
    f = open(file, 'r')
    lines = f.readlines()

    table = dict()
    counter = 0
    for line in lines:
        line_to_list = line.split('\t')
        id = line_to_list[0]
        count = int(line_to_list[1])
        tfidf = int(line_to_list[2])
        name = line_to_list[3]
        mass_string = line_to_list[4]
        if mass_string == 'NaN':
            mass = 777
        else:
            mass = float(mass_string)
        logP = float(line_to_list[5])
        superterms = line_to_list[6]
        table[id] = {"Count": count, "TFIDF": tfidf, "Name": name, "Mass": mass, "logP": logP, "Superterms": superterms}

    return table

def read_file(file):
    '''
    This function reads a file and adds the information + id in a dictionary. This is added to another dictionary as value, and the term that describes the
    information (e.g. "Names") as key.
    '''
    f = open(file, 'r')
    lines = f.readlines()
    id_to_name = dict()
    for line in lines:
        line_to_list = line.split('\t')
        id = line_to_list[0]
        name = line_to_list[1].strip()
        id_to_name[id] = name

    return id_to_name

def create_array(table, normalization):
    '''
    This function recieves the table dictionary and returns numpy arrays with the mass, logP and names of the chemicals.
    These values and names are repeated corresponding to the count value.
    '''
    x = []
    y = []
    ids = []
    counts_total = 0
    tfidf_total = 0
    for id in table.keys():
        counts_total += table[id]['Count']
        tfidf_total += table[id]['TFIDF']

    normalizing_factor = counts_total / tfidf_total

    list_counts = []
    for id in table.keys():
        if normalization:
            count = table[id]["TFIDF"]
            count = int(count * normalizing_factor)+1
            list_counts.append(count)
        else:
            count = table[id]["Count"]
        mass = table[id]["Mass"]
        logP = table[id]["logP"]
        x.extend([logP]*count)
        y.extend([mass]* count)
        ids.extend([id]*count)

    return np.asarray(x), np.asarray(y), ids

def hexbin(x, y, ids, size, aspect_scale, orientation):
    '''
    This function recieves x, y coordinates arrays and turns these into q, r hexagon coordinates returned in a pandas dataframe.
    Additionally, the center x, y coordinates for every hexagon are added.
    '''
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=aspect_scale)
    df = pd.DataFrame(dict(r=r, q=q))
    df['ids'] = ids
    df = df.groupby(['q', 'r'])['ids'].apply(list).reset_index(name='ids')
    counts = []
    for row in df['ids']:
        counts.append(len(row))

    df['counts'] = counts
    return df

def find_most_common(list_names, tooltip_count):
    '''
    This function recieves a list of (repeated) chemical ids. It returns a list of the tooltip_count=X most common chemical ids and their counts.
    In some hexagons there are less than 5 chemicals, if so empty names/counts are returned.
    '''
    uniques = list(set(list_names))
    name_to_count = {}

    # count the amount of appearances for every name in this hexagon
    for name in uniques:
        count = list_names.count(name)
        name_to_count[name] = count

    # sort by count
    from collections import Counter
    k = Counter(name_to_count)

    # find 5 highest values
    high = k.most_common(tooltip_count)

    # make a list of X names, if there are no more names in high, append empty strings as names ...
    # so that something can be added to the columns when this list is returned
    most_common_names = []
    for i in range(tooltip_count):
        if (i+1) <= len(high):
            name_and_count = high[i]
        else:
            name_and_count = ('', '')
        most_common_names.append(name_and_count)

    return most_common_names

def add_counts(df, table):
    '''
    This function recieves a pandas dataframe with q, r coordinates and names as columns.
    For every row (hexagon), the total amount of chemicals is counted and the 5 most common chemicals are put in seperate columns.
    '''
    # Amount of ids shown in the tooltip (e.g. with TOOLTIP_COUNT=3 the 3 most common chebi ids will be shown)
    # Currently, max value is 3 (see html code for hover tooltip, probably needs to be changed as well when changing this value)
    TOOLTIP_COUNT = 3

    total_counts = []

    # make columns
    ids = [[] for _ in range(TOOLTIP_COUNT)]
    counts = [[] for _ in range(TOOLTIP_COUNT)]
    names = [[] for _ in range(TOOLTIP_COUNT)]

    for row in df['ids']:
        # total amount of hits
        total_counts.append(len(row))

        # count 5 most common unique hits (names), returns list
        most_common_chemicals = find_most_common(row, TOOLTIP_COUNT)

        # add names to their corresponding columns (most common name in first list, 2nd most common name in second list etc.)
        for i in range(TOOLTIP_COUNT):
            id = most_common_chemicals[i][0]
            count = most_common_chemicals[i][1]
            try:
                name = table[id]["Name"]
            except:
                name = ""

            ids[i].append(id)
            counts[i].append(count)
            names[i].append(name)

    #add lists to df as columns
    for i in range(TOOLTIP_COUNT):
        column_name = "id"+str(i+1)
        df[column_name] = ids[i]
        column_name = "count"+str(i+1)
        df[column_name] = counts[i]
        column_name = "name"+str(i+1)
        df[column_name] = names[i]

    df['counts'] = total_counts # original count values
    df['scaling'] = total_counts # column used for root scaling

    return df

def get_blur(x,y,sigma_x,sigma_y):
    '''
    This function recieves x, y values and sigma x, sigma y values and returns the calculated blur value.
    See https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    '''

    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def add_gaussian_blur(df, blur_max, step_size):
    '''
    Function:
    This function adds gaussian blur to the plot by going through all hexagons and applying the kernel to the neighbouring hexagons.
    To speed up the process, the pandas dataframe is put in a python dictionary, so to quickly find the neighbouring hexagon coordinates using the coordinates as keys.
    After bluring, the dictionary is put in a pandas dataframe again, and this dataframe is returned.

    Columns:
    sd_x values for calculating blur values are used to name a column with counts that result from blurring with that specific sd_x.
    This makes selecting correct sd_x column easy with the slider code, because the slider returns values that represent sd_x values.
    ColumnDataSource does not accept intergers as column names, so sd_x column names are changed to string.
    sd_x columns contain lists, in a list the first value is the normal counts, the second value is the tfidf count.
    '''
    # columns for blur
    blur_columns_count = int(blur_max / step_size)

    # add blur columns to dataframe
    blur_list = []

    # prepare column names for the new dataframe
    column_names = []

    for column_name in df.columns:
        column_names.append(column_name)

    for i in range(0, blur_columns_count+1):
        sd_x = i*step_size
        if sd_x%1 == 0:
            sd_x = int(sd_x) # this to get the names of the column names equal to the bokeh slider values, e.g. slider returns 0 or 1 when python returns 0.0 or 1.0
        blur_list.append(sd_x)
        column_name = str(sd_x)
        column_names.append(column_name)

    # construct kernel with distances
    kernel_distances = {(-5,2):(7.5, sqrt(3)/2),(-5,3):(7.5, sqrt(3)/2),
    (-4,1):(6, sqrt(3)),(-4,2):(6, 0),(-4,3):(6,sqrt(3)),
    (-3,0):(4.5, 3*sqrt(3)/2),(-3,1):(4.5, sqrt(3)/2),(-3,2):(4.5, sqrt(3)/2),(-3,3):(4.5, 3*sqrt(3)/2),
    (-2,-1):(3, 2*sqrt(3)),(-2,0):(3, sqrt(3)),(-2,1):(3, 0),(-2,2):(3, sqrt(3)),(-2,3):(3, 2*sqrt(3)),
    (-1,-1):(1.5, 3*sqrt(3)/2),(-1,0):(1.5, sqrt(3)/2),(-1,1):(1.5, sqrt(3)/2),(-1,2):(1.5, 3*sqrt(3)/2),
    (0,-2):(0, 2*sqrt(3)),(0,-1):(0, sqrt(3)),(0,1):(0, sqrt(3)),(0,2):(0, 2*sqrt(3)),
    (1,-2):(1.5, 3*sqrt(3)/2),(1,-1):(1.5, sqrt(3)/2),(1,0):(1.5, sqrt(3)/2),(1,1):(1.5, 3*sqrt(3)/2),(2,-3):(3, 2*sqrt(3)),
    (2,-2):(3, sqrt(3)),(2,-1):(3, 0),(2,0):(3, sqrt(3)),(2,1):(3, 2*sqrt(3)),
    (3,-3):(4.5, 3*sqrt(3)/2),(3,-2):(4.5, sqrt(3)/2),(3,-1):(4.5, sqrt(3)/2),(3,0):(4.5, 3*sqrt(3)/2),
    (4,-3):(6, sqrt(3)),(4,-2):(6, 0),(4,-1):(6, sqrt(3)),
    (5,-3):(7.5, sqrt(3)/2),(5,-2):(7.5, sqrt(3)/2)}

    # make kernel for blur
    kernel_blur = {}
    for sd_x in blur_list:
        sd_y = sd_x / 2

        # calculate blur for all the keys in the kernel, replace distance with blur factor
        for key in kernel_distances:
            distance = kernel_distances [key]
            x = distance[0]
            y = distance[1]
            if sd_x == 0:
                blur = 0
            else:
                blur = get_blur(x, y, sd_x, sd_y)
            try:
                kernel_blur[key][sd_x] = blur
            except:
                kernel_blur[key] = {sd_x: blur}

    # Remake df in python dictionary
    coords_to_count = dict() # keeps all original coordinates and count values
    coords_to_blur = dict() # will have new coordinates and new count (blur) values

    # go through all coordinates in df
    for index, row in df.iterrows():
        q = row['q']
        r = row['r']
        counts = row['counts']
        tfidf = row['tfidf']
        coords_to_blur[(q,r)] = {}
        coords_to_count[(q,r)] = {'counts': counts, 'tfidf': tfidf}

        # add every column from the dataframe to the dictionary
        for i in range(len(df.columns)):
            column_name = df.columns[i]
            column_value = row[column_name]
            coords_to_blur[(q,r)][column_name] = column_value

    # put original count values in sd_x column
    for sd_x in blur_list:
        column_name = str(sd_x)
        # Add blur
        for coord in coords_to_count.keys():
            q = coord[0]
            r = coord[1]
            counts = coords_to_count[coord]['counts']
            tfidf = coords_to_count[coord]['tfidf']
            coords_to_blur[(q, r)][column_name] = [counts, tfidf]


    # do blur for every blur_column
    for sd_x in blur_list:
        column_name = str(sd_x)
        # Add blur
        for coord in coords_to_count.keys():
            q = coord[0]
            r = coord[1]
            counts = coords_to_count[coord]['counts']
            tfidf = coords_to_count[coord]['tfidf']

            # For every coordinate, apply kernel (to surrounding hexagons)
            for key in kernel_blur.keys():
                q_new = q + key[0]
                r_new = r + key[1]
                factor = kernel_blur[key][sd_x]
                counts_blur = counts*factor
                tfidf_blur = tfidf*factor

                # Check if surrounding hexagon exists and apply blur
                try:
                    coords_to_blur[(q_new, r_new)][column_name][0] += counts_blur
                    coords_to_blur[(q_new, r_new)][column_name][1] += tfidf_blur
                except:
                    # If hexagon does not exist, make a new hexagon
                    # All columns (see keys) need to be added here or else Bokeh will crash, but most values will be left empty
                    if sd_x == 0:
                        # make the row
                        coords_to_blur[(q_new, r_new)] = {'q': q_new, 'r': r_new, column_name: [counts_blur, tfidf_blur], 'counts': 0, 'tfidf': 0, 'scaling': 0, 'ids': [],
                        'id1': '', 'count1': '', 'name1': '',
                        'id2': '', 'count2': '', 'name2': '',
                        'id3': '', 'count3': '', 'name3': '',
                        'id1_tfidf': '', 'count1_tfidf': '', 'name1_tfidf': '',
                        'id2_tfidf': '', 'count2_tfidf': '', 'name2_tfidf': '',
                        'id3_tfidf': '', 'count3_tfidf': '', 'name3_tfidf': ''}
                    else:
                        # row should already be there, only need to add column value
                        coords_to_blur[(q_new, r_new)][column_name] = [counts_blur, tfidf_blur]
                        # coords_to_blur[(q_new, r_new)][column_name] = [counts_blur, tfidf_blur]


    df_gaussian = pd.DataFrame()

    # return dictionary to pandas dataframe
    for column_name in column_names:
        # add all dictionary keys as columns (lists) to the dataframe
        column_list = []
        for coords in coords_to_blur.keys():
            column_value = coords_to_blur[coords][column_name]
            column_list.append(column_value)
        # add column to the df
        df_gaussian[column_name] = column_list

    df_gaussian = df_gaussian.drop("ids", axis=1)
    return df_gaussian

def merge_dataframes(df_normal, df_tfidf):
    df_normal['tfidf'] = df_tfidf['counts']
    df_normal['id1_tfidf'] = df_tfidf['id1']
    df_normal['id2_tfidf'] = df_tfidf['id2']
    df_normal['id3_tfidf'] = df_tfidf['id3']
    df_normal['count1_tfidf'] = df_tfidf['count1']
    df_normal['count2_tfidf'] = df_tfidf['count2']
    df_normal['count3_tfidf'] = df_tfidf['count3']
    df_normal['name1_tfidf'] = df_tfidf['name1']
    df_normal['name2_tfidf'] = df_tfidf['name2']
    df_normal['name3_tfidf'] = df_tfidf['name3']
    return df_normal

def create_class_source(table, term, size, ratio, orientation, superterm):
    id_to_name = read_file('files/ChEBI2Names.tsv')
    x = []
    y = []
    counter = 0

    for id in table.keys():
        superterms = table[id]["Superterms"]
        superterms = superterms.strip()
        superterms = superterms.strip("[]")
        superterms = superterms.split(',')
        for class_id in superterms:
            class_id = class_id.replace("\"", '')
            class_id = class_id.replace("\'", '')
            class_id = class_id.strip()
            try:
                name = id_to_name[class_id]
            except:
                print([class_id])

            if name == superterm:
                mass = table[id]["Mass"]
                logP = table[id]["logP"]
                x.extend([logP])
                y.extend([mass])
    x = np.asarray(x)
    y = np.asarray(y)
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)
    df = pd.DataFrame(dict(r=r, q=q))
    source = ColumnDataSource(df)

    return source

def create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE):
    # create array
    x, y, ids = create_array(table, normalization=False)
    xn, yn, idsn = create_array(table, normalization=True)

    # title
    title = 'Hexbin plot for '+str(len(x))+' annotated chemicals with query '+str(term)

    # make dataframe with counts
    df = hexbin(x, y, ids, size, aspect_scale=ratio, orientation=orientation)
    dfn = hexbin(xn, yn, idsn, size, aspect_scale=ratio, orientation=orientation)

    # add counts for chemicals for hexagon
    df = add_counts(df, table)
    dfn = add_counts(dfn, table)
    df = merge_dataframes(df, dfn)
    df = add_gaussian_blur(df, BLUR_MAX, BLUR_STEP_SIZE)
    source = ColumnDataSource(df)

    return source, title

def make_plot_sources(tables, size, ratio, orientation, blur_max, step_size):
    '''
    This function is used to make lists of datasources for the plot that can each be selected through use of the widgets.
    These data sources are sources with or without normalization, and a range of applied blur (each blur is on a new dataset).

    This funcion recieves:
    - tables: results from the query search (chebi ids, mass, logp etc.), these will be used to create the x, y coordinates
    - size, ratio, orientation: plot properties used to make the hexagonal bins
    - blur_max, step_size: blur properties used to create gaussian blur

    This function returns a dictionary with terms (names of the tables) as keys and lists of data sources as values.

    '''

    term_to_source = dict()
    term_to_metadata = dict()
    options = []

    for term in tables.keys():
        print('making dataframe for %s ...' % term)
        # options for multi select
        options.append((term, term))
        metadata = tables[term]['metadata']
        metadata = return_html(metadata)
        table = tables[term]['table']
        # tables[term]['table'] = [] # delete from memory ## DOENST WORK


        #create array
        x, y, ids = create_array(table, normalization=False)
        xn, yn, idsn = create_array(table, normalization=True)

        # title
        title = 'Hexbin plot for '+str(len(x))+' annotated chemicals with query '+str(term)


        df = hexbin(x, y, ids, size, aspect_scale=ratio, orientation=orientation)
        dfn = hexbin(xn, yn, idsn, size, aspect_scale=ratio, orientation=orientation)

        # add counts for chemicals for hexagon
        df = add_counts(df, table)
        dfn = add_counts(dfn, table)

        # MERGE DFS

        df = add_gaussian_blur(df, blur_max, step_size)
        source = ColumnDataSource(df)

        # Add litst to dictionary
        term_to_source[term] = {'source': source, 'title': title}
        term_to_metadata[term] = metadata


    return term_to_source, term_to_metadata, options

def return_JS_code(widget):
    '''
    This function recieves a string that indicates the widget.
    It returns a JavaScript callback code corresponding to the widget.
    '''
    if widget == 'tooltips':
        code = """
            <style>
            table, td, th {
              border-collapse: collapse;
              border: 1px solid #dddddd;
              padding: 2px;
              table-layout: fixed;
              height: 20px;
            }

            tr:nth-child(even) {
              background-color: #dddddd;
            }
            </style>

            <table>
              <col width="100">
              <col width="80">
              <col width="220">
              <tr>
                <th>Total Counts</th>
                <th>@counts</th>
                <th>($x, $y)</th>
              </tr>
              <tr style="color: #fff; background: black;">
                <th>Chebi ID</th>
                <th>Count</th>
                <th>Name</th>
              </tr>
              <tr>
                <th>@id1</th>
                <th>@count1</th>
                <th>@name1</th>
              </tr>
              <tr>
                <th>@id2</th>
                <th>@count2</th>
                <th>@name2</th>
              </tr>
              <tr>
                <th>@id3</th>
                <th>@count3</th>
                <th>@name3</th>
              </tr>
            </table>

            """
    elif widget == 'tooltips_tfidf':
        code = """
            <style>
            table, td, th {
              border-collapse: collapse;
              border: 1px solid #dddddd;
              padding: 2px;
              table-layout: fixed;
              height: 20px;
            }

            tr:nth-child(even) {
              background-color: #dddddd;
            }
            </style>

            <table>
              <col width="100">
              <col width="80">
              <col width="220">
              <tr>
                <th>Total Counts</th>
                <th>@tfidf</th>
                <th>($x, $y)</th>
              </tr>
              <tr style="color: #fff; background: black;">
                <th>Chebi ID</th>
                <th>Count</th>
                <th>Name</th>
              </tr>
              <tr>
                <th>@id1_tfidf</th>
                <th>@count1_tfidf</th>
                <th>@name1_tfidf</th>
              </tr>
              <tr>
                <th>@id2_tfidf</th>
                <th>@count2_tfidf</th>
                <th>@name2_tfidf</th>
              </tr>
              <tr>
                <th>@id3_tfidf</th>
                <th>@count3_tfidf</th>
                <th>@name3_tfidf</th>
              </tr>
            </table>

            """

    elif widget == 'slider1':
        code = """
            var mapper = mapper
            var source_data = source.data;
            var f = cb_obj.value
            var sd_x =  Math.round(slider2.value)
            var sd_x = String(sd_x)
            var checkbox = checkbox

            if (checkbox.active.lenght == 1) {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][1], 1/f)
                }
            } else {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][0], 1/f)
                }
            }

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['scaling'])

            // apply changes
            source.change.emit();
            """

    elif widget == 'slider2':
        code = """
            var source_data = source.data;
            var sd_x =  cb_obj.value;
            var sd_x = String(sd_x)
            var f = slider1.value;
            var checkbox = checkbox;
            var mapper = mapper;

            // apply scaling

            if (checkbox.active.length == 1) {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][1], 1/f)
                }
            } else {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][0], 1/f)
                }
            }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['scaling'])

            // apply changes
            source.change.emit();
            """

    elif widget == 'rbg':
        code = """
            var active = cb_obj.active;
            var p = p;
            var Viridis256 = Viridis256;
            var Greys256 = Greys256;
            var class_hex = class_hex;
            var term_to_class = term_to_class;
            var term = multi_select.value[0]

            if (active == 1){
            mapper.transform.palette = Greys256
            p.background_fill_color = '#000000'
            if (term_to_class[term]['show_class']){class_hex.glyph.fill_color = 'red'}
            }

            if (active == 0){
            mapper.transform.palette = Viridis256
            p.background_fill_color = '#440154'
            if(term_to_class[term]['show_class']){class_hex.glyph.fill_color = 'pink'}
            }

            """

    elif widget == 'checkbox':
        code = """
            var hover = hover;
            var tooltips = tooltips
            var tooltips_tfidf = tooltips_tfidf
            var source_data = source.data;
            var active = cb_obj.active
            var f = slider1.value
            var sd_x =  Math.round(slider2.value)
            var sd_x = String(sd_x)
            var mapper = mapper

            // apply scaling

            if (active.length == 1) {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][1], 1/f)
                }

            hover.tooltips = tooltips_tfidf

            } else {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][0], 1/f)
                }

            hover.tooltips = tooltips
            }

            // max value for mapper
            mapper.transform.high = Math.max.apply(Math, source_data['scaling'])

            source.change.emit();
            """

    elif widget == 'hover':
        code = """
            var tooltips = document.getElementsByClassName("bk-tooltip");
            for (var i = 0, len = tooltips.length; i < len; i ++) {
                tooltips[i].style.top = ""; // unset what bokeh.js sets
                tooltips[i].style.left = "";
                tooltips[i].style.bottom = "150px";
                tooltips[i].style.left = "575px";
                tooltips[i].style.width = "500px";
            }
            """

    elif widget == 'button':
        code = """
            var term_to_metadata = term_to_metadata
            var term = multi_select.value[0]
            var metadata = term_to_metadata[term]
            var wnd = window.open("about:blank", "", "_blank");
            wnd.document.write(metadata)
            """

    elif widget == 'multi_select':
        code = """
            var source_data = source.data;
            var term_to_source = term_to_source;
            var term = cb_obj.value[0];
            var f = slider1.value;
            var sd_x =  Math.round(slider2.value)
            var sd_x = String(sd_x)
            var checkbox = checkbox;
            var p = p;
            var mapper = mapper;
            var metadata = metadata;
            var class_hex = class_hex;
            var checkbox_class = checkbox_class
            var term_to_class = term_to_class;

            // select new_data
            var new_data = term_to_source[term]['source'].data

            // replace current data
            for (key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // select correct blur column, apply saturation
            if (checkbox.active.length == 1) {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][1], 1/f)
                }
            } else {

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data[sd_x][i][0], 1/f)
                }
            }

            // new title
            var title = term_to_source[term]['title']
            p.title.text = title

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['scaling'])

            // class
            if (term_to_class[term]['show_class']){
                class_hex.visible = false
                class_hex.source = term_to_class[term]['source']
                checkbox_class.active = []
            }
            source.change.emit();
            """
    elif widget == 'class':
        code = """
            var term_to_class = term_to_class;
            var multi_select = multi_select;
            var class_hex = class_hex;
            var active = cb_obj.active;

            if (active.length == 1) {
                class_hex.visible = true
            } else {
                class_hex.visible = false
            }


            """
    return code

def return_html(metadata):
    '''
    This function recieves the metadata. The metadata is put in the html string and this string returned.
    '''
    html_content = """
    <HTML>
    <HEAD>
    <TITLE>%s</TITLE>
    </HEAD>
    <BODY BGCOLOR="FFFFFF">
    %s
    <br>%s
    <br>%s
    <br>%s
    <br>%s
    <br>%s
    </BODY>
    </HTML>
    """ % (metadata[0], metadata[0], metadata[1], metadata[2], metadata[3], metadata[4], metadata[5])
    return html_content

def plot(tables, output_filename, xmin, xmax, ymin, ymax, superterm):
    '''
    This is the plot function that uses Bokeh functions and widgets to make an interactive hexagon plot.

    This function recieves:
    - tables: dictionary with tables used to create arrays of repeated x, y coordinates (depending on the counts) for the hexagon plot.
    - output_filename: filename of .html output in the plots folder

    The coordinate arrays are used to create a pandas dataframe with Bokeh functions. This dataframe contains the q, r coordinates and counts used to plot the
    hexagons. To this dataframe, extra information is added (e.g. most common chemicals), which is displayed in the hover tooltip.

    Gaussian blur is added to copies of this dataframe and given as input to the Bokeh slider widget.
    Other widgets are added as well, for saturation, normalization etc. Bokeh allows to customize these widges with javascript code.

    The hexagon plot is saved as a .html file and also shown in the browser.
    '''

    file_name = 'plots/'+str(output_filename)+'.html'
    output_file(file_name)

    # Blur and saturation values
    BLUR_MAX = 4
    BLUR_STEP_SIZE = 0.25
    SATURATION_MAX = 5
    SATURATION_STEP_SIZE = 0.25

    # Hexagon plot properties
    SIZE_HEXAGONS = 10
    orientation = 'flattop' #bokeh alows 2 different hexagon orientations which also influences hexagon size calculations, but we currently have only calculated blur distances for this orientation
    ratio = ((ymax-ymin) / (xmax-xmin) )
    size = SIZE_HEXAGONS / ratio
    hexagon_height = sqrt(3) * size
    hexagon_height = hexagon_height*ratio

    # make figure
    p = figure(x_range = [xmin, xmax],y_range=[ymin-(hexagon_height/2), ymax],
               tools="wheel_zoom,reset,save", background_fill_color= '#440154')

    p.grid.visible = False
    p.xaxis.axis_label = "log(P)"
    p.yaxis.axis_label = "mass in Da"
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'

    # term_to_source, term_to_metadata, options = make_plot_sources(tables, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)

    # source for widgets
    term_to_source = dict()
    term_to_class = dict()
    term_to_metadata = dict()
    options = []

    for term in tables.keys():
        options.append((term, term))
        table = tables[term]['table']
        if superterm:
            source = create_class_source(table, term, size, ratio, orientation, superterm)
            term_to_class[term] = {}
            term_to_class[term]['show_class'] = True
            term_to_class[term]['source'] = source
        else:
            term_to_class[term] = {'show_class': False}
        source, title = create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
        metadata = return_html(tables[term]['metadata'])
        term_to_source[term] = {'source': source, 'title': title}
        term_to_metadata[term] = metadata

    # hex = p.hex_tile(q='q', r="r", size=size, line_color=None, source=source, aspect_scale=ratio,orientation=orientation,
    # fill_color='pink' )

    # show(p)

    # make default souce for plot, this is the first source shown in the plot, and also works like a container. Old data is thrown out and new data is thrown in.
    default_term = list(tables.keys())[0] # pick the first one
    metadata = tables[default_term]['metadata']
    metadata = return_html(metadata)
    table = tables[default_term]['table']
    source, title = create_data_source(table, default_term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
    p.title.text = title

    # color mapper
    mapper = linear_cmap('scaling', 'Viridis256', 0, max(source.data['scaling']))

    # plot
    hex = p.hex_tile(q="q", r="r", size=size, line_color=None, source=source, aspect_scale=ratio, orientation=orientation,
           fill_color=mapper)

    if superterm:
        source_class = term_to_class[default_term]['source']
        class_hex = p.hex_tile(q='q', r="r", size=size, line_color=None, source=source_class, aspect_scale=ratio,orientation=orientation,
            fill_color='pink', fill_alpha=0.7)
        class_hex.visible = False

    # HOVER
    TOOLTIPS = return_JS_code('tooltips')
    TOOLTIPS_tfidf = return_JS_code('tooltips_tfidf')
    code_callback_hover = return_JS_code('hover')
    callback_hover = CustomJS(code=code_callback_hover)
    hover = HoverTool(renderers=[hex], tooltips=TOOLTIPS, callback=callback_hover, show_arrow=False)
    p.add_tools(hover)

    # WIDGETS
    slider1 = Slider(start=1, end=SATURATION_MAX, value=1, step=SATURATION_STEP_SIZE, title="Saturation", width=100)
    slider2 = Slider(start=0, end=BLUR_MAX, value=0, step=BLUR_STEP_SIZE, title="Blur", width=100)
    checkbox = CheckboxGroup(labels=["TFIDF"], active=[])
    radio_button_group = RadioGroup(labels=["Viridis256", "Greys256"], active=0)
    button = Button(label="Metadata",button_type="default", width=100)
    multi_select = MultiSelect(title=output_filename, value=[default_term],options=options, width=100, height=300)
    if superterm:
        label = "Show "+str(superterm)
        checkbox_class = CheckboxGroup(labels=[label], active=[])

    # WIDGETS CODE FOR CALLBACK
    code_callback_slider1 = return_JS_code('slider1')
    code_callback_slider2 = return_JS_code('slider2')
    code_callback_checkbox = return_JS_code('checkbox')
    code_callback_rbg = return_JS_code('rbg')
    code_callback_button = return_JS_code('button')
    code_callback_ms = return_JS_code('multi_select')
    if superterm:
        code_callback_class = return_JS_code('class')

    # WIDGETS CALLBACK
    callback_slider1 = CustomJS(args={'source': source, 'mapper': mapper, 'slider2': slider2, 'checkbox': checkbox}, code=code_callback_slider1)
    callback_slider2 = CustomJS(args={'source': source, 'mapper': mapper, 'slider1': slider1, 'checkbox': checkbox}, code=code_callback_slider2)
    callback_checkbox = CustomJS(args={'source': source, 'slider1': slider1, 'slider2': slider2, 'mapper': mapper, 'hover': hover, 'tooltips': TOOLTIPS, 'tooltips_tfidf': TOOLTIPS_tfidf}, code=code_callback_checkbox)
    callback_radio_button_group = CustomJS(args={'p': p, 'multi_select': multi_select, 'mapper': mapper, 'term_to_class': term_to_class, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
    callback_button = CustomJS(args={'term_to_metadata': term_to_metadata, 'multi_select': multi_select},code=code_callback_button)
    callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source,  'term_to_class': term_to_class, 'checkbox': checkbox, 'slider2': slider2, 'slider1': slider1, 'p': p, 'mapper': mapper}, code=code_callback_ms)
    if superterm:
        callback_radio_button_group = CustomJS(args={'p': p, 'multi_select': multi_select, 'class_hex': class_hex, 'term_to_class': term_to_class, 'mapper': mapper, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
        callback_class = CustomJS(args={'multi_select': multi_select, 'term_to_class': term_to_class, 'class_hex': class_hex}, code=code_callback_class)
        callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source, 'checkbox': checkbox, 'slider2': slider2, 'slider1': slider1, 'p': p, 'mapper': mapper,
        'checkbox_class': checkbox_class, 'class_hex': class_hex, 'term_to_class': term_to_class}, code=code_callback_ms)

    # WIDGETS INTERACTION
    slider1.js_on_change('value', callback_slider1)
    slider2.js_on_change('value', callback_slider2)
    checkbox.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio_button_group)
    button.js_on_event(events.ButtonClick, callback_button)
    multi_select.js_on_change("value", callback_ms)
    if superterm:
        checkbox_class.js_on_change('active', callback_class)

    # LAYOUT
    if superterm:
        layout = row(multi_select, p, column(slider1, slider2, checkbox, checkbox_class, radio_button_group, button))
    else:
        layout = row(multi_select, p, column(slider1, slider2, checkbox, radio_button_group, button))

    show(layout)


def get_tables(files):
    '''
    This function recieves a list of files, imports these files and their corresponding metadata, and saves them in a dictionary.
    '''
    tables = dict()
    for file in files:
        table = import_table(file)
        term = file.split('\\')[1].split('_')[0]
        metadata_file = 'metadata/'+str(term)+'.txt'
        metadata_lines = open(metadata_file, 'r')
        metadata = metadata_lines.readlines()

        tables[term] = {'table': table, 'metadata': metadata}
    return tables

def get_files(folder):
    '''
    This function recieves the input folder and returns a list of file paths in this folder.
    '''
    files = []
    for file in listdir(folder):
        path = str(folder)+'\\'+str(file)
        files.append(path)
    return files

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input_folder', dest='input_folder', help='[i] to select input folder')
    parser.add_argument('-o', required=True, metavar='output_filename', dest='output_filename', help='[o] to name .html output filename (a \'plots\' folder is required!)')
    parser.add_argument('-xmin', required=False, metavar='xmin', dest='xmin', help='[xmin] to select x axis minimum (logP), default is -5')
    parser.add_argument('-xmax', required=False, metavar='xmax', dest='xmax', help='[xmax] to select x axis maximum (logP), default is 10')
    parser.add_argument('-ymax', required=False, metavar='ymax', dest='ymax', help='[ymax] to select y axix maximum (mass in Da), default is 1600')
    parser.add_argument('-class', required=False, metavar='superterm', dest='superterm', help='[c] to select a class to be shown in the plot with a click on the class button')
    arguments = parser.parse_args()
    return arguments

def main():
    startTime = datetime.now()
    args = parser()
    folder = args.input_folder
    output_filename = args.output_filename
    xmin = args.xmin
    xmax = args.xmax
    ymax = args.ymax

    # default settings
    if not xmin:
        xmin = -5
    else:
        xmin = float(args.xmin)
    if not xmax:
        xmax = 10
    else:
        xmax = float(args.xmax)
    if not ymax:
        ymax = 1600
    else:
        ymax = float(args.ymax)
    ymin = 0 # cannot be negative

    files = get_files(folder)
    tables = get_tables(files)

    #
    plot(tables, output_filename, xmin, xmax, ymin, ymax, args.superterm)
    # print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
