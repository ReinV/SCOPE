#!/usr/bin/python

import numpy as np
import pandas as pd
import argparse
from math import exp
from math import sqrt
from datetime import datetime

from bokeh import events
from bokeh.io import output_file, show
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
        mass = float(line_to_list[4])
        logP = float(line_to_list[5])

        table[id] = {"Count": count, "TFIDF": tfidf, "Name": name, "Mass": mass, "logP": logP}

    return table

def create_array(table, normalisation):
    '''
    This function recieves the table dictionary and returns numpy arrays with the mass, logP and names of the chemicals.
    These values and names are repeated corresponding to the count value.
    '''
    x = []
    y = []
    ids = []
    for id in table.keys():
        if normalisation:
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

def find_most_common(list_names):
    '''
    This function recieves a list of (repeated) chemical ids. It returns a list of the 5 most common chemical ids and their counts.
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
    high = k.most_common(5)

    # EXPLAIN TUPPLE?
    # make a list of 5 names, if there are no more names in high, append empty strings as names ...
    # so that something can be added to the columns when this list is returned
    most_common_names = []
    for i in range(5):
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
    total_counts = []

    # make columns
    ids = [[] for _ in range(5)]
    counts = [[] for _ in range(5)]
    names = [[] for _ in range(5)]

    for row in df['ids']:
        # total amount of hits
        total_counts.append(len(row))

        # count 5 most common unique hits (names), returns list
        most_common_chemicals = find_most_common(row)

        # add names to their corresponding columns (most common name in first list, 2nd most common name in second list etc.)
        for i in range(5):
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
    for i in range(5):
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

def construct_source_list(df):
    '''
    This function recieves a dataframe and adds gaussian blur to copies of this dataframe for a range of blur values.
    It returns a list with these dataframes ordered from no blur to high blur.
    '''
    source_list = [ColumnDataSource(data=df)]
    for i in range(1, 17):
        sd_x = i / 4
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

    # calculate blur for all the keys in the kernel, replace distance with blur factor
    for key in kernel:
        distance = kernel[key]
        x = distance[0]
        y = distance[1]
        blur = get_blur(x, y, sd_x, sd_y)
        kernel[key] = blur

    # make normal python dict
    coords_to_count = dict()
    coords_to_blur = dict()

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

    # go through all coordinates
    for coord in coords_to_count.keys():
        q = coord[0]
        r = coord[1]
        counts = coords_to_count[coord]
        # for every coordinate, apply kernel (to surrounding hexagons)
        for key in kernel.keys():
            q_new = q + key[0]
            r_new = r + key[1]
            factor = kernel[key]
            # check if surrounding hexagon exists and apply blur
            try:
                coords_to_blur[(q_new, r_new)]['blur'] += counts*factor
                coords_to_blur[(q_new, r_new)]['scaling'] += counts*factor
            except:
                # if hexagon does not exist, make a new hexagon (all columns (keys) need to be added, but most values will be empty)
                blur_new = counts*factor
                coords_to_blur[(q_new, r_new)] = {'q': q_new, 'r': r_new, 'x': '-', 'y': '-','counts': 0, 'blur': blur_new, 'scaling': blur_new, 'ids': [],
                'id1': '', 'count1': '', 'name1': '',
                'id2': '', 'count2': '', 'name2': '',
                'id3': '', 'count3': '', 'name3': '',
                'id4': '', 'count4': '', 'name4': '',
                'id5': '', 'count5': '', 'name5': ''}

    # new dataframe with gaussian blur
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
      <tr>
        <th>@id4</th>
        <th>@count4</th>
        <th>@name4</th>
      </tr>
      <tr>
        <th>@id5</th>
        <th>@count5</th>
        <th>@name5</th>
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
            var value = value*4;
            var mapper = mapper;
            var power = slider1.value;
            var source_list = source_list;
            var source_norm_list = source_norm_list;
            var source_data = source.data;
            var checkbox = checkbox

            if (checkbox.value == 0){
            var new_data = source__norm_list[value].data;
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
            var active = cb_obj.active
            var source_data = source.data;
            var source_norm_list = source_norm_list;
            var source_list = source_list;
            var value1 = slider1.value
            var value2 = slider2.value
            var mapper = mapper



            if (active.length == 1 ) {

            // replace source data with tfidf data
            new_data = source_norm_list[value2*4].data
            for (key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                    }
                    }
            } else {

            // replace source data with data without tfidf
            new_data = source_list[value2*4].data
            for (key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                    }
                    }
            }

            // apply current saturation level
            for (var i = 0; i < source_data['blur'].length; i++) {
                source_data['scaling'][i] = Math.pow(source_data['blur'][i], 1/value1)
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
            var html_content = html_content
            var wnd = window.open("about:blank", "", "_blank")
            wnd.document.write(html_content)
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

def plot(term, table, metadata):
    '''
    This is the plot function that uses Bokeh functions and widgets to make an interactive hexagon plot.

    This function recieves:
    - term: string used for the title of the plot and naming the output file.
    - table: dictionary used to create arrays of repeated x, y coordinates (depending on the counts) for the hexagon plot.
    - metadata: string of info that will be shown with the plot

    The coordinate arrays are used to create a pandas dataframe with Bokeh functions. This dataframe contains the q, r coordinates and counts used to plot the
    hexagons. To this dataframe, extra information is added (e.g. most common chemicals), which is displayed in the hover tooltip.

    Gaussian blur is added to copies of this dataframe and given as input to the Bokeh slider widget.
    Other widgets are added as well, for saturation, normalisation etc. Bokeh allows to customize these widges with javascript code.

    The hexagon plot is saved as a .html file and also shown in the browser.
    '''

    file_name = 'plots/'+str(term)+'.html'
    output_file(file_name)

    x, y, ids = create_array(table, normalisation=False)
    xn, yn, idsn = create_array(table, normalisation=True)

    # hexagon plot properties
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

    # make dataframe with counts
    df = hexbin(x, y, ids, size, aspect_scale=ratio, orientation=orientation)
    dfn = hexbin(xn, yn, idsn, size, aspect_scale=ratio, orientation=orientation)

    # add counts for hexagon
    df = add_counts(df, table)
    dfn = add_counts(dfn, table)

    # add gausian blur, add columns with names/counts, put dfs with blur in list
    source_list = construct_source_list(df)
    source_norm_list = construct_source_list(dfn)

    # source for plot
    source = ColumnDataSource(data=df)

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
    slider1 = Slider(start=1, end=5, value=1, step=0.25, title="Saturation", width=100)
    slider2 = Slider(start=0, end=4, value=0, step=0.25, title="Blur", width=100)
    checkbox = CheckboxGroup(labels=["TFIDF"], active=[])
    radio_button_group = RadioGroup(labels=["Viridis256", "Greys256"], active=0)
    button = Button(label="Metadata",button_type="default", width=100)

    # WIDGETS CODE FOR CALLBACK
    code_callback_slider1 = return_code('slider1')
    code_callback_slider2 = return_code('slider2')
    code_callback_checkbox = return_code('checkbox')
    code_callback_rbg = return_code('rbg')
    code_callback_button = return_code('button')

    # HTML
    html_content = return_html(metadata)

    # WIDGETS CALLBACK
    callback_slider1 = CustomJS(args={'source': source, 'mapper': mapper}, code=code_callback_slider1)
    callback_slider2 = CustomJS(args={'source': source, 'source_list': source_list, 'source_norm_list': source_norm_list,
        'slider1': slider1, 'checkbox': checkbox, 'mapper': mapper}, code=code_callback_slider2)
    callback_checkbox = CustomJS(args={'source': source, 'source_norm_list': source_norm_list, 'source_list': source_list, 'mapper': mapper, 'slider1': slider1, 'slider2': slider2},
        code=code_callback_checkbox)
    callback_radio_button_group = CustomJS(args={'p': p, 'mapper': mapper, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
    button_callback = CustomJS(args={'html_content': html_content},code=code_callback_button)

    # WIDGETS INTERACTION
    slider1.js_on_change('value', callback_slider1)
    slider2.js_on_change('value', callback_slider2)
    checkbox.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio_button_group)
    button.js_on_event(events.ButtonClick, button_callback)

    # LAYOUT
    layout = row(p, column(slider1, slider2, checkbox, radio_button_group, button))

    show(layout)


def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    arguments = parser.parse_args()
    return arguments

def main():
    startTime = datetime.now()
    args = parser()
    file = args.input_file
    table = import_table(file)
    term = file.split('\\')[1].split('_')[0]

    metadata_file = 'metadata/'+str(term)+'.txt'
    metadata_lines = open(metadata_file, 'r')
    metadata = metadata_lines.readlines()

    plot(term, table, metadata)

    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
