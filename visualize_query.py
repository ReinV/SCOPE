#!/usr/bin/python

import numpy as np
import pandas as pd
import argparse
from math import exp
from math import sqrt
from datetime import datetime
from os import listdir
import sys

from bokeh import events
from bokeh.io import output_file, show, save
from bokeh.models import HoverTool
from bokeh.models import CustomJS
from bokeh.models import HoverTool
from bokeh.models import CustomJS, ColumnDataSource, Slider, CheckboxGroup, RadioGroup, Button
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

        table[id] = {"Count": count, "TFIDF": tfidf, "Name": name, "Mass": mass, "logP": logP}

    return table

def create_array(table, normalization):
    '''
    This function recieves the table dictionary and returns numpy arrays with the mass, logP and names of the chemicals.
    These values and names are repeated corresponding to the count value.
    '''
    x = []
    y = []
    ids = []
    for id in table.keys():
        if normalization:
            count = table[id]["TFIDF"]
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
    x, y = axial_to_cartesian(q, r, size, orientation=orientation, aspect_scale=aspect_scale)
    df = pd.DataFrame(dict(r=r, q=q, x=x, y=y))
    df['ids'] = ids
    df = df.groupby(['q', 'r', 'x', 'y'])['ids'].apply(list).reset_index(name='ids')

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
    # Currently, max value is 3 (see html code for hover tooltip)
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
    df['blur'] = total_counts # blur values
    df['scaling'] = total_counts # column used for root scaling

    return df

def construct_source_list(df, blur_max, step_size):
    '''
    This function recieves a dataframe and adds gaussian blur to copies of this dataframe for a range of blur values.
    It returns a list with these dataframes ordered from no blur to high blur.
    '''
    source_list = [ColumnDataSource(df)]
    end = int(blur_max/step_size)+1
    for i in range(1, end):
        sd_x = i*step_size
        sd_y = sd_x / 2
        df_new = add_gaussian_blur(df, sd_x, sd_y)
        source = ColumnDataSource(data=df_new)
        source_list.append(source)

    return source_list

def get_blur(x,y,sigma_x,sigma_y):
    '''
    This function recieves x, y values and sigma x, sigma y values and returns the calculated blur value.
    See https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    '''
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def add_gaussian_blur(df, sd_x, sd_y):
    '''
    This function adds gaussian blur to the plot by going through all hexagons and applying the kernel to the neighbouring hexagons.
    To speed up the process, the panda dataframe is put in a python dictionary, so to quickly find the neighbouring hexagon coordinates using the coordinates as keys.
    After bluring, the dictionary is put in a panda dataframe, and this dataframe is returned.
    '''
    # construct kernel with distances
    kernel = {(-5,2):(7.5, sqrt(3)/2),(-5,3):(7.5, sqrt(3)/2),
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


    # make normal python dict
    coords_to_count = dict()
    coords_to_blur = dict()

    # calculate blur for all the keys in the kernel, replace distance with blur factor
    for key in kernel:
        distance = kernel[key]
        x = distance[0]
        y = distance[1]
        blur = get_blur(x, y, sd_x, sd_y)
        kernel[key] = blur

    # go through all coordinates
    for index, row in df.iterrows():
        q = row['q']
        r = row['r']
        counts = row['counts']
        coords_to_blur[(q,r)] = {}

        # add dataframe to a python dictionary
        coords_to_count[(q,r)] = counts
        for i in range(len(df.columns)):
            column_name = df.columns[i]
            column_value = row[i]
            coords_to_blur[(q,r)][column_name] = column_value

    # Add blur
    for coord in coords_to_count.keys():
        q = coord[0]
        r = coord[1]
        counts = coords_to_count[coord]
        # For every coordinate, apply kernel (to surrounding hexagons)
        for key in kernel.keys():
            q_new = q + key[0]
            r_new = r + key[1]
            factor = kernel[key]
            # Check if surrounding hexagon exists and apply blur
            try:
                coords_to_blur[(q_new, r_new)]['blur'] += counts*factor
                coords_to_blur[(q_new, r_new)]['scaling'] += counts*factor
            except:
                # If hexagon does not exist, make a new hexagon
                # All columns (keys) need to be added here or else Bokeh will crash, but most values will be left empty
                blur_new = counts*factor
                coords_to_blur[(q_new, r_new)] = {'q': q_new, 'r': r_new, 'x': '-', 'y': '-','counts': 0, 'blur': blur_new, 'scaling': blur_new, 'ids': [],
                'id1': '', 'count1': '', 'name1': '',
                'id2': '', 'count2': '', 'name2': '',
                'id3': '', 'count3': '', 'name3': ''}

    df_gaussian = pd.DataFrame()

    # add all dictionary keys as columns (lists) to the dataframe
    for column in df.columns:
        column_list = []
        for coords in coords_to_blur.keys():
            column_value = coords_to_blur[coords][column]
            column_list.append(column_value)
        # add column to the df
        df_gaussian[column] = column_list

    return df_gaussian

def return_tooltip():
    '''
    This function returns the HTML code for the hover tooltip.
    '''
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
        <th>(@x,@y)</th>
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
    return code

def return_code(widget):
    '''
    This function recieves a string that indicates the widget.
    It returns a javascript callback code corresponding to the widget.
    '''
    if widget == 'slider1':
        code = """
            var mapper = mapper
            var data = source.data;
            var f = cb_obj.value
            var c = data['blur']
            var s = data['scaling']

            // apply scaling
            for (var i = 0; i < c.length; i++) {
                s[i] = Math.pow(c[i], 1/f)
            }

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, s)

            // apply changes
            source.change.emit();
            """

    elif widget == 'slider2':
        code = """
            var value = cb_obj.value;
            var step_size = step_size
            var value = Math.round(value/step_size)
            var mapper = mapper;
            var power = slider1.value
            var source_data = source.data;
            var source_list = source_list;
            var source_norm_list = source_norm_list;
            var checkbox = checkbox;

            if (checkbox.active.length == 1){
                var new_data = source_norm_list[value].data;
            } else {
                var new_data = source_list[value].data;
            }

            // apply blur
            for (key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // apply current saturation level
            for (var i = 0; i < source_data['blur'].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data['blur'][i], 1/power)
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

            if (active == 1){
            mapper.transform.palette = Greys256
            p.background_fill_color = '#000000'
            }

            if (active == 0){
            mapper.transform.palette = Viridis256
            p.background_fill_color = '#440154'
            }

            """

    elif widget == 'checkbox':
        code = """

            var source_data = source.data;
            var source_list = source_list;
            var source_norm_list = source_norm_list;
            var active = cb_obj.active
            var power = slider1.value
            var value = slider2.value
            var step_size = step_size
            var value = Math.round(value/step_size)
            var mapper = mapper

            if (active.length == 1 ) {
                    var new_data = source_norm_list[value].data
                } else {
                    var new_data = source_list[value].data
                }

            // replace dataset
            for (key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // apply current saturation level
            for (var i = 0; i < source_data['blur'].length; i++) {
            source_data['scaling'][i] = Math.pow(source_data['blur'][i], 1/power)
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
            var metadata = metadata;
            var wnd = window.open("about:blank", "", "_blank");
            wnd.document.write(metadata)
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

def make_plot_source(table, term, size, ratio, orientation, blur_max, step_size):
    '''
    This function is used to make lists of datasources for the plot that can each be selected through use of the widgets.
    These data sources are sources with or without normalization, and a range of applied blur (each blur is on a new dataset).

    This funcion recieves:
    - tables: results from the query search (chebi ids, mass, logp etc.), these will be used to create the x, y coordinates
    - size, ratio, orientation: plot properties used to make the hexagonal bins
    - blur_max, step_size: blur properties used to create gaussian blur

    This function returns a dictionary with terms (names of the tables) as keys and lists of data sources as values.

    '''

    #create array
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

    # add blur
    source = construct_source_list(df, blur_max, step_size)
    source_norm = construct_source_list(dfn, blur_max, step_size)

    return source, source_norm, title

def plot(table, metadata, term):
    '''
    This is the plot function that uses Bokeh functions and widgets to make an interactive hexagon plot.

    This function recieves:
    - table: dictionary with chebi identifiers as keys and their (normalized) counts, mass, logP etc. as values
    - metadata: metadata (date, query etc.)
    - term: term used to name the output file

    The coordinate arrays are used to create a pandas dataframe with Bokeh functions. This dataframe contains the q, r coordinates and counts used to plot the
    hexagons. To this dataframe, extra information is added (e.g. most common chemicals), which is displayed in the hover tooltip.

    Gaussian blur is added to copies of this dataframe and given as input to the Bokeh slider widget.
    Other widgets are added as well, for saturation, normalization etc. Bokeh allows to customize these widges with javascript code.

    The hexagon plot is saved as a .html file and also shown in the browser.
    '''

    file_name = 'plots/'+str(term)+'.html'
    output_file(file_name)

    # Blur and saturation values
    BLUR_MAX = 3
    BLUR_STEP_SIZE = 0.5
    SATURATION_MAX = 5
    SATURATION_STEP_SIZE = 0.25

    # First, create array for plot properties ( ratio, size of hexagons etc.)
    x, y, ids = create_array(table, normalization=False)

    # Hexagon plot properties
    length = len(x)
    orientation = 'flattop'
    ratio = ( (max(y)-min(y)) / (max(x)-min(x)) )
    size = 10 / ratio
    h = sqrt(3) * size
    h = h*ratio
    title = 'Hexbin plot for '+str(length)+' annotated chemicals with query '+str(term)

    # make figure
    p = figure(title=title, x_range = [min(x)-0.5, max(x)+0.5],y_range=[0-(h/2),max(y)+100],
               tools="wheel_zoom,reset,save", background_fill_color= '#440154')

    p.grid.visible = False
    p.xaxis.axis_label = "log(P)"
    p.yaxis.axis_label = "mass in Da"
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'

    # source for plot with blur for dataframes with and without normalization
    source_list, source_norm_list, title = make_plot_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)

    # start source for plot, this is the source that is first displayed in the hexagon figure and is replaced whenever changes are made with the widgets
    df = hexbin(x, y, ids, size, aspect_scale=ratio, orientation=orientation)
    df = add_counts(df, table)
    source = ColumnDataSource(df)

    # color mapper
    mapper = linear_cmap('scaling', 'Viridis256', 0, max(source.data['scaling']))

    # plot
    hex = p.hex_tile(q="q", r="r", size=size, line_color=None, source=source, aspect_scale=ratio, orientation=orientation,
           fill_color=mapper)

    # HOVER
    TOOLTIPS = return_tooltip()
    code_callback_hover = return_code('hover')
    callback_hover = CustomJS(code=code_callback_hover)
    hover = HoverTool(tooltips=TOOLTIPS, callback=callback_hover, show_arrow=False)
    p.add_tools(hover)

    # WIDGETS
    slider1 = Slider(start=1, end=SATURATION_MAX, value=1, step=SATURATION_STEP_SIZE, title="Saturation", width=100)
    slider2 = Slider(start=0, end=BLUR_MAX, value=0, step=BLUR_STEP_SIZE, title="Blur", width=100)
    checkbox = CheckboxGroup(labels=["TFIDF"], active=[])
    radio_button_group = RadioGroup(labels=["Viridis256", "Greys256"], active=0)
    button = Button(label="Metadata",button_type="default", width=100)

    # WIDGETS CODE FOR CALLBACK
    code_callback_slider1 = return_code('slider1')
    code_callback_slider2 = return_code('slider2')
    code_callback_checkbox = return_code('checkbox')
    code_callback_rbg = return_code('rbg')
    code_callback_button = return_code('button')

    # WIDGETS CALLBACK
    callback_slider1 = CustomJS(args={'source': source, 'mapper': mapper}, code=code_callback_slider1)
    callback_slider2 = CustomJS(args={'source': source, 'mapper': mapper, 'slider1': slider1, 'checkbox': checkbox, 'source_list': source_list,
    'source_norm_list': source_norm_list, 'step_size': BLUR_STEP_SIZE}, code=code_callback_slider2)
    callback_checkbox = CustomJS(args={'source': source, 'source_list': source_list, 'source_norm_list': source_norm_list, 'step_size': BLUR_STEP_SIZE,
    'slider1': slider1, 'slider2': slider2, 'mapper': mapper}, code=code_callback_checkbox)
    callback_radio_button_group = CustomJS(args={'p': p, 'mapper': mapper, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
    callback_button = CustomJS(args={'metadata': metadata},code=code_callback_button)

    # WIDGETS INTERACTION
    slider1.js_on_change('value', callback_slider1)
    slider2.js_on_change('value', callback_slider2)
    checkbox.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio_button_group)
    button.js_on_event(events.ButtonClick, callback_button)

    # LAYOUT
    layout = row(p, column(slider1, slider2, checkbox, radio_button_group, button))

    save(layout)

def get_tables(files):
    '''
    This function recieves a list of files, imports these files and their corresponding metadata, and saves them in a dictionary.
    '''
    tables = dict()
    for file in files:
        table = import_table(file)
        term = file.split('\\')[1].split('_')[0]
        metadata = get_metadata(term)
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

def get_metadata(term):
    metadata_file = 'metadata/'+str(term)+'.txt'
    metadata_lines = open(metadata_file, 'r')
    metadata = metadata_lines.readlines()
    return metadata

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input', dest='input', help='[i] to select input')
    parser.add_argument('-t', required=True, metavar='input_type', dest='input_type', help='[t] to select input type: folder or file')
    arguments = parser.parse_args()
    return arguments

def main():
    startTime = datetime.now()
    args = parser()
    input = args.input
    input_type = args.input_type

    if input_type == 'folder':
        files = get_files(input)
        tables = get_tables(files)
        # plot all files in the folder
        for term in tables.keys():
            output_filename = term
            table = table[term]
            metadata = table[metadata]
            print('plotting for %s ...' % term)
            plot(table, metadata, term)

    elif input_type == 'file':
        # plot the file
        table = import_table(input)
        term = input.split('\\')[1].split('_')[0]
        metadata = get_metadata(term)
        plot(table, metadata, term)

    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
