#!/usr/bin/python

import numpy as np
import pandas as pd
import argparse

from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import hexbin
from bokeh.util.hex import cartesian_to_axial
from math import exp
from math import sqrt

from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, Slider

from datetime import datetime

def hex_to_bin(x, y, size, aspect_scale, orientation):
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=aspect_scale)
    df = pd.DataFrame(dict(r=r, q=q))
    return df

def import_table(file):
    f = open(file, 'r')
    lines = f.readlines()

    table = dict()
    for line in lines:
        line_to_list = line.split('\t')
        id = line_to_list[0]
        count = int(line_to_list[1])
        name = line_to_list[2]
        mass = float(line_to_list[3])
        logP = float(line_to_list[4])

        table[id] = {"Count": count, "Name": name, "Mass": mass, "logP": logP}

    return table

def create_array(table):
    x = []
    y = []
    names = []
    for id in table.keys():
        mass = table[id]["Mass"]
        logP = table[id]["logP"]
        count = table[id]["Count"]
        name = table[id]["Name"]
        x.extend([logP]*count)
        y.extend([mass]* count)
        names.extend([name]*count)
    return np.asarray(x), np.asarray(y), names

def find_most_common(list_names):
    uniques = list(set(list_names))
    most_common_name = ""
    count = 0
    for name in uniques:
        new_count = list_names.count(name)
        if new_count > count:
            count = new_count
            most_common_name = name
    return most_common_name

def transform_df(df):
    counts = []
    most_common_list = []
    for element in df['names']:
        counts.append(len(element))
        name = find_most_common(element)
        most_common_list.append(name)
    df['names'] = most_common_list
    df['counts'] = counts
    # column that will be used for color scaling
    # these numbers are different when gaussian blur is applyied but we want to keep the original counts column for the hover
    df['color_scaling'] = counts
    return df

def plot(x, y, names, term):

    from bokeh.models import CustomJS
    from bokeh.models import HoverTool


    callback_hover = CustomJS(code="""
        var tooltips = document.getElementsByClassName("bk-tooltip");
        for (var i = 0, len = tooltips.length; i < len; i ++) {
            tooltips[i].style.top = ""; // unset what bokeh.js sets
            tooltips[i].style.left = "";
            tooltips[i].style.bottom = "-60px";
            tooltips[i].style.left = "270px";
        }
        """)

    length = len(x)
    orientation = 'flattop'
    size = 10
    ratio = ((max(y)-min(y)) / (max(x)-min(x)) )

    if orientation == 'flattop':
        size = size / ratio

    title = 'Hexbin plot for '+str(length)+' annotated chemicals with query '+str(term)

    # make figure
    p = figure(title=title, x_range = [min(x)-1, max(x)],y_range=[min(y)-100,max(y)],
               tools="wheel_zoom,reset", background_fill_color= '#440154')
    p.grid.visible = False
    p.xaxis.axis_label = "log(P)"
    p.yaxis.axis_label = "mass in Da"
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'

    # make dataframe
    df = hex_to_bin(x, y, size, ratio, orientation)


    # add names to dataframe column 'names'
    df['names'] = names

    # for every bin (row with r, q coordinates), group all names in a list
    df = df.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')

    # for every bin, count the total amount of chemicals and return most common checmial
    df_00 = transform_df(df)

    print('applying gaussian blur...')
    # gaussian blur is added with sd_x and sd_y as arguments
    df_05 = add_gaussian_blur(df_00, 0.5, 0.25)
    df_10 = add_gaussian_blur(df_00, 1, 0.5)
    df_15 = add_gaussian_blur(df_00, 1.5, 0.75)
    df_20 = add_gaussian_blur(df_00, 2, 1)
    df_25 = add_gaussian_blur(df_00, 2.5, 1.25)
    df_30 = add_gaussian_blur(df_00, 3, 1.5)

    # for slider
    source1 = ColumnDataSource(data=df_00)
    source2 = ColumnDataSource(data=df_05)
    source3 = ColumnDataSource(data=df_10)
    source4 = ColumnDataSource(data=df_15)
    source5 = ColumnDataSource(data=df_20)
    source6 = ColumnDataSource(data=df_25)
    source7 = ColumnDataSource(data=df_30)
    source = ColumnDataSource(data=df_00)

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=source,aspect_scale=ratio, orientation=orientation,
           fill_color=linear_cmap('color_scaling', 'Viridis256', 0, max(source.data['color_scaling'])))

    # tooltips = TOOLTIPS
    # [("count", "@counts"),("name", "@names")]
    hover = HoverTool(tooltips=[("q","@q"),("r","@r"),("count", "@counts"),("name", "@names")],callback=callback_hover)
    #  mode="mouse", point_policy="follow_mouse", renderers=[r2]
    p.add_tools(hover)
    #
    callback = CustomJS(args={'source': source, 'source1': source1, 'source2': source2, 'source3': source3, 'source4': source4, 'source5': source5, 'source6': source6, 'source7': source7},
    code="""
        var f = cb_obj.value;
        var sdata = source.data;
        var data1 = source1.data;
        var data2 = source2.data;
        var data3 = source3.data;
        var data4 = source4.data;
        var data5 = source5.data;
        var data6 = source6.data;
        var data7 = source7.data;

        console.log(data2);

        for (key in data1) {console.log(key);}

        if (f == "0") {
        for (key in data1) {
            sdata[key] = [];
            for (i=0;i<data1[key].length;i++){
            sdata[key].push(data1[key][i]);
            }
        }
        } else if (f == "0.5") {
        for (key in data2) {
            sdata[key] = [];
            for (i=0;i<data2[key].length;i++){
            sdata[key].push(data2[key][i]);
            }
        }
        } else if (f == "1") {
        for (key in data3) {
            sdata[key] = [];
            for (i=0;i<data3[key].length;i++){
            sdata[key].push(data3[key][i]);
            }
        }
        } else if (f == "1.5") {
        for (key in data4) {
            sdata[key] = [];
            for (i=0;i<data4[key].length;i++){
            sdata[key].push(data4[key][i]);
            }
        }
        } else if (f == "2") {
        for (key in data5) {
            sdata[key] = [];
            for (i=0;i<data5[key].length;i++){
            sdata[key].push(data5[key][i]);
            }
        }
        } else if (f == "2.5") {
        for (key in data6) {
            sdata[key] = [];
            for (i=0;i<data6[key].length;i++){
            sdata[key].push(data6[key][i]);
            }
        }
        } else if (f == "3") {
        for (key in data7) {
            sdata[key] = [];
            for (i=0;i<data7[key].length;i++){
            sdata[key].push(data7[key][i]);
            }
        }
        }

        source.change.emit();
        """)

    # Add the callback to the select widget.
    # This executes each time the selected option changes


    slider = Slider(start=0, end=3, value=0, step=0.5, title="sigma_x")
    slider.js_on_change('value', callback)

    layout = column(slider, p)

    show(layout)


def get_blur(x,y,sigma_x,sigma_y):
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def add_gaussian_blur(df, sd_x, sd_y):
    '''
    This function adds gaussian blur to the plot by going through all hexagons and applying the kernel to the neighbouring hexagons.
    '''

    # new def
    df_gaussian = df

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

    # loop through rows
    for index, row in df_gaussian.iterrows():
        q = row['q']
        r = row['r']
        names = row['names']
        counts = row['color_scaling']

        # for every q,r coordinate, apply kernel to neighbouring hexagons
        for key in kernel.keys():
            q_new = q + key[0]
            r_new = r + key[1]
            blur = kernel[key]

            # check if neighbouring hexagon exists
            try:
                # look up counts for new coordinates
                counts_old = df_gaussian[(df_gaussian.q==q_new) & (df_gaussian.r==r_new)]['color_scaling'].iloc[0]
                # get index for this hexagon
                index_new = df_gaussian.index[(df_gaussian.q==q_new) & (df_gaussian.r==r_new)]
                index_new = int(index_new[0])
                # add blur to the color_scaling counts
                df_gaussian.at[index_new, 'color_scaling'] = counts_old + counts*blur

            except:
                # if hexagon does not exist, append it to the dataframe with new coordinates, no counts, but with added blur
                df_gaussian = df_gaussian.append({'q': q_new, 'r': r_new, 'counts': 0, 'color_scaling': (blur * counts)}, ignore_index = True)

    return df_gaussian

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

    x, y, names = create_array(table)
    plot(x, y, names, term)

    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
