#!/usr/bin/python

import numpy as np
import pandas as pd
import argparse

from bokeh.models import CustomJS, CheckboxGroup, RadioButtonGroup
from bokeh.layouts import column
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import hexbin
from bokeh.util.hex import cartesian_to_axial
from math import exp
from math import sqrt

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
        tfidf = int(line_to_list[2])
        name = line_to_list[3]
        mass = float(line_to_list[4])
        logP = float(line_to_list[5])

        table[id] = {"Count": count, "tfidf": tfidf, "Name": name, "Mass": mass, "logP": logP}

    return table

def create_array(table, normalization):
    x = []
    y = []
    names = []
    for id in table.keys():
        if normalization:
            count = table[id]["tfidf"]
        else:
            count = table[id]["Count"]
        mass = table[id]["Mass"]
        logP = table[id]["logP"]
        name = table[id]["Name"]
        x.extend([logP]*count)
        y.extend([mass]* count)
        names.extend([name]*count)
    return np.asarray(x), np.asarray(y), names

def find_most_common(list_names):
    uniques = list(set(list_names))
    name_to_count = {}

    # count the amount of appearances for every name in this hexagon,
    for name in uniques:
        count = list_names.count(name)
        name_to_count[name] = count

    # sort by count
    from collections import Counter
    k = Counter(name_to_count)

    # find 5 highest values
    high = k.most_common(5)

    # make a list of 5 names, if there are no more names in high, append empty strings as names ..
    # so that something can be added to the columns when this list is returned
    most_common_names = []
    for i in range(5):
        if (i+1) <= len(high):
            name_and_count = high[i]
        else:
            name_and_count = ''
        most_common_names.append(name_and_count)

    return most_common_names

def transform_df(df):
    counts = []
    # make columns
    columns = [[] for _ in range(5)]
    for row in df['names']:
        # total amount of hits
        counts.append(len(row))

        # count 5 most common unique hits (names), returns list
        most_common_names = find_most_common(row)

        # add names to their corresponding columns (most common name in first list, 2nd most common name in second list etc.)
        for i in range(5):
            name = most_common_names[i]
            columns[i].append(name)

    #add lists to df as columns
    for i in range(5):
        column_name = "names"+str(i+1)
        df[column_name] = columns[i]
    print(df.columns)
    df['counts'] = counts

    return df

def plot(table, term):

    #create coordinates
    x, y, names = create_array(table,normalization=False)
    xn, yn, names_n = create_array(table,normalization=True)
    # hexagon properties
    orientation = 'flattop'
    size = 10
    ratio = ((max(y)-min(y)) / (max(x)-min(x)) )
    if orientation == 'flattop':
        size = size / ratio

    # make hexagon dataframe
    df = hex_to_bin(x, y, size, ratio, orientation)
    df_normalised = hex_to_bin(xn, yn, size, ratio, orientation)
    # add names to dataframe column 'names'
    df['names'] = names
    df_normalised['names'] = names_n
    # for every bin (row with r, q coordinates), group all names in a list
    df = df.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    df_normalised = df_normalised.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    # for every bin, count the total amount of chemicals and return most common checmial
    df = transform_df(df)
    df_normalised = transform_df(df_normalised)
    # add small gaussian blur with sdx = 2 and sdy = 1
    # df = add_gaussian_blur(df)
    # df_normalised = add_gaussian_blur(df_normalised)
    # make the df in column data source for the callback stuff later in the function
    # one extra copy of source wihout normalization for the callback of the checkbox group
    source = ColumnDataSource(data=df)
    source_std = ColumnDataSource(data=df)
    source_norm = ColumnDataSource(data=df_normalised)
    # mappers for scale of colors
    mapper1 = linear_cmap('counts', 'Viridis256', 0, max(source.data['counts']))
    mapper2 = log_cmap('counts', 'Viridis256', 0, max(source.data['counts']))
    # mapper3 = linear_cmap('gaussian_blur', 'Viridis256', 0, max(source.data['gaussian_blur']))
    # mapper4 = log_cmap('gaussian_blur', 'Viridis256', 0, max(source.data['gaussian_blur']))

    #plot title
    length = len(x)
    title = 'Hexbin plot for '+str(length)+' annotated chemicals with query '+str(term)
    #(snap_to_data, follow_mouse, none) ) poitn policy
    # make figure
    p = figure(title=title, x_range = [min(x)-1, max(x)],y_range=[min(y)-100,max(y)],
    tools="wheel_zoom,reset", background_fill_color= '#440154')
    p.grid.visible = False
    p.xaxis.axis_label = "log(P)"
    p.yaxis.axis_label = "mass in Da"
    p.xaxis.axis_label_text_font_style = 'normal'
    p.yaxis.axis_label_text_font_style = 'normal'
    # make hextile plot
    hex = p.hex_tile(q="q", r="r", size=size, line_color=None, source=source,aspect_scale=ratio, orientation=orientation,
               hover_color="pink", hover_alpha=0.8, fill_color=mapper1, legend='Viridis256')
    #hover position
    callback_hover = CustomJS(code="""
        var tooltips = document.getElementsByClassName("bk-tooltip");
        for (var i = 0, len = tooltips.length; i < len; i ++) {
            tooltips[i].style.top = ""; // unset what bokeh.js sets
            tooltips[i].style.left = "";
            tooltips[i].style.bottom = "280px";
            tooltips[i].style.left = "575px";
            tooltips[i].style.width = "300px";
        }
        """)
    hover = HoverTool(tooltips=[("total count", "@counts"),("#1", "@names1"),("#2", "@names2"),("#3", "@names3"),("#4", "@names4"),("#5", "@names5")],callback=callback_hover)
    p.add_tools(hover)

    #widgets

    # radio button for log or linear (replace mapper)
    radio_button_group = RadioButtonGroup(
        labels=["Linear scale", "Logarithmic scale"], active=0)

    # checkbox for TFIDF (replace dataframe) and gaussian blur (replace mapper)
    checkbox_group = CheckboxGroup(
        labels=["TF IDF", "Gaussian blur"], active=[])

    # define callbacks for the widgest
    # , 'mapper3': mapper3, 'mapper4': mapper4
    callback_radio = CustomJS(args={"hex": hex, 'checkbox': checkbox_group, 'mapper1': mapper1, 'mapper2': mapper2}, code="""
        var active = cb_obj.active;
        var checkbox = checkbox.active

        // check if "Logarithmic" is selected
        if (active == 1){

        // check if "gaussian blur" is selected with checkbox
        if (checkbox.includes(1)){
        hex.glyph.fill_color = mapper4
        } else {
        hex.glyph.fill_color = mapper2
        }
        } else {
        if (checkbox.includes(1)){
        hex.glyph.fill_color = mapper3
        } else{
        hex.glyph.fill_color = mapper1
        }
        }
        hex.change.emit()
        """)
# , "mapper3": mapper3, "mapper4": mapper4
    callback_checkbox = CustomJS(args={"source": source, "source_std": source_std, "source_norm": source_norm, "hex": hex, "radio": radio_button_group,
    "mapper1": mapper1, "mapper2": mapper2}, code="""
        var active = cb_obj.active;
        var data_source = source.data;
        var data_std = source_std.data;
        var data_norm = source_norm.data;
        var hex = hex;
        var radio = radio;

        // check if option "TF IDF" is checked for normalization
        // if so, replace normal data source with normalized data source (or else the reverse)

        if (active.includes(0)){
            for (key in data_norm) {
                data_source[key] = [];
                for (i=0;i<data_norm[key].length;i++){
                data_source[key].push(data_norm[key][i]);
        }}
        } else {
            for (key in data_std) {
                data_source[key] = [];
                for (i=0;i<data_std[key].length;i++){
                data_source[key].push(data_std[key][i]);
        }}}

        // check if option "gaussian blur" is checked
        // if so, scale hex plot color according to gaussian blur in stead of normal counts
        // also check the current scale (active=0: linear, or active==1:log) from the radio button group

        if (active.includes(1)){
            if (radio.active == 0){
            hex.glyph.fill_color = mapper3

            } else {
            hex.glyph.fill_color = mapper4
        }
        } else {
            if (radio.active == 0){
                hex.glyph.fill_color = mapper1
            } else{
                hex.glyph.fill_color = mapper2
            }
        }

        hex.change.emit()
        source.change.emit()
        """)

    checkbox_group.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio)

    layout = column(checkbox_group, radio_button_group, p)

    show(layout)



def get_blur(x,y):
    sigma_x = 2
    sigma_y = 1
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def add_gaussian_blur(df):
    '''
    This function adds gaussian blur to the plot by checking every possible coordinate with a gaussian kernel and calculating a new value.
    The kernel has a center, inside ring and outside ring with attached weights.
    '''

    # new def
    df_gaussian = df
    df_gaussian['gaussian_blur'] = df_gaussian['counts']

    # construct kernel
    # kernel = {(0,-1):get_blur(0,sqrt(3)),(1,-1):get_blur(1.5,sqrt(3)/2),(1,0): get_blur(1.5,sqrt(3)/2),(2,-1): get_blur(3,0),
    # (0,1):get_blur(0,sqrt(3)),(-1,1) :get_blur(1.5,sqrt(3)/2),(-1,0): get_blur(1.5,sqrt(3)/2),(-2,1):get_blur(3,0)}

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
        blur = get_blur(x, y)
        kernel[key] = blur

    # loop through rows
    for index, row in df_gaussian.iterrows():
        q = row['q']
        r = row['r']
        names = row['names']
        counts = row['gaussian_blur']

        for key in kernel.keys():
            q_new = q + key[0]
            r_new = r + key[1]
            blur = kernel[key]
            try:
                counts_old = df_gaussian[(df_gaussian.q==q_new) & (df_gaussian.r==r_new)]['gaussian_blur'].iloc[0]
                index_new = df_gaussian.index[(df_gaussian.q==q_new) & (df_gaussian.r==r_new)]
                index_new = int(index_new[0])
                df_gaussian.at[index_new, 'gaussian_blur'] = counts_old + counts*blur
            except:
                df_gaussian = df_gaussian.append({'q': q_new, 'r': r_new, 'counts': 0, 'gaussian_blur': (blur * counts)}, ignore_index = True)

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
    plot(table, term)

    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
