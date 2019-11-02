#!/usr/bin/python

import numpy as np
import pandas as pd
import math
import argparse

from collections import Counter
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
from bokeh.palettes import Reds9 as red
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import hexbin
from bokeh.util.hex import cartesian_to_axial


def hex_to_bin(x, y, size, aspect_scale, orientation="pointytop"):
    # pd = import_required('pandas','hexbin requires pandas to be installed')

    q, r = cartesian_to_axial(x, y, size, orientation, aspect_scale=aspect_scale)
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
        superterms = line_to_list[5]
        superterms = superterms.strip()
        superterms = superterms.strip("\"")
        superterms = superterms.strip('[]')
        superterms = superterms.replace('\'','')
        superterms = superterms.split(', ')

        table[id] = {"Count": count, "Name": name, "Mass": mass, "logP": logP, "Superterms": superterms}

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
    return df

def most_abundant_classes(table):
    classes = dict()
    for key in table.keys():
        for superterm in table[key]["Superterms"]:
            superterm = superterm.strip("\"")
            try:
                classes[superterm] += 1
            except:
                classes[superterm] = 1

    k = Counter(classes)

    # Finding 10 highest values
    high = k.most_common(100)

    return high

def select_data(table, classes):
    data = dict()
    for id in table.keys():
        superterms = table[id]["Superterms"]
        if classes in superterms:
            data[id] = table[id]
        else:
            pass

    return data

def plot_class(x, y, names, xs, ys, names_s):

    red.reverse()
    white8 = 8*['#FFFFFF']
    grey8 = 8*['#808080']

    ratio = (max(y)-min(y)) / (max(x)-min(x))
    size = 20
    # (max(y)-min(y)) / (max(x)-min(x))

    # print(ratio)
    p = figure(title="Hexbin for many points", match_aspect=True, aspect_scale = ratio, x_range = [min(x), max(x)],y_range=[min(y),max(y)],
               tools="wheel_zoom,reset", background_fill_color= '#FFFFFF')
    p.grid.visible = False

    alpha = 0.8
    # r, bins = p.hexbin(x, y, size=0.5, hover_color="pink", hover_alpha=0.8)

    df = hex_to_bin(x, y, size, ratio)
    dfs = hex_to_bin(xs, ys, size, ratio)
    # df_m2 = hex_to_bin(a, b, size=0.2)

    df['names'] = names
    dfs['names'] = names_s
    # df_m2['names'] = names?

    df = df.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    dfs = dfs.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')

    df = transform_df(df)
    dfs = transform_df(dfs)

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df,aspect_scale=ratio,
               fill_color=log_cmap('counts', grey8, 0, max(df.counts)))

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=dfs, aspect_scale=ratio,
                fill_color=log_cmap('counts', red, 0, max(dfs.counts)))

    # r, bins = p.hexbin(x, y, size=0.5)

    hover = HoverTool(tooltips=[("count", "@counts"),("name", "@names")])
    # renderers=[r2]
    p.add_tools(hover)
    #
    output_file("plots/RPC_class_red_test1.html")

    show(p)
    return

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=False, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    parser.add_argument('-c', required=False, metavar='classes', dest='classes', help='[c] to select class to show in plot')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    file = args.input_file
    classes = args.classes
    table = import_table(file)

    selected_table = select_data(table, classes)#  select only class from the data
    x, y, names = create_array(table)
    xs, ys, names_s = create_array(selected_table)
    plot_class(x, y, names, xs, ys, names_s)


    # high = most_abundant_classes(table)
    # # print(high[200:210])
    # x, y, names = create_array_no_tfidf(table)

if __name__ == '__main__':
    main()
