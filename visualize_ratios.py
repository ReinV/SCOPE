#!/usr/bin/python

import numpy as np
import pandas as pd
import math
import argparse

from collections import Counter
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
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
        tfidf = int(line_to_list[2])
        name = line_to_list[3]
        mass = float(line_to_list[4])
        logP = float(line_to_list[5])
        superterms = line_to_list[6]
        superterms = superterms.strip()
        superterms = superterms.strip("\"")
        superterms = superterms.strip('[]')
        superterms = superterms.replace('\'','')
        superterms = superterms.split(', ')

        table[id] = {"Count": count,"tfidf": tfidf, "Name": name, "Mass": mass, "logP": logP, "Superterms": superterms}

    return table

def create_array(table):
    x = []
    y = []
    names = []
    for id in table.keys():
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

def calculate_ratios(df, dfs, lb):

    df_ratio_low = pd.DataFrame()
    df_ratio_high = pd.DataFrame()

    q_start = min(min(df['q']),min(dfs['q']))
    q_end = max(max(df['q']),max(dfs['q']))
    r_start = min(min(df['r']),min(dfs['r']))
    r_end = max(max(df['r']),max(dfs['r']))

    q_list_low = []
    r_list_low = []
    ratio_list_low = []

    q_list_high = []
    r_list_high = []
    ratio_list_high = []

    for q in range(q_start, q_end+1):
        for r in range(r_start, r_end+1):
            try:
                count = df[(df.q==q) & (df.r==r)]['counts'].iloc[0] # check if count exists for these coordinates
                name = df[(df.q==q) & (df.r==r)]['names'].iloc[0]
            except:
                count = 0
                name = ""
            try:
                count_s = dfs[(dfs.q==q) & (dfs.r==r)]['counts'].iloc[0]
                name_s = df[(df.q==q) & (df.r==r)]['names'].iloc[0]
            except:
                count_s = 0
                name_s = ""

            if count == 0 and count_s == 0:
                pass
            else:
                if count < lb:
                    count = lb
                if count_s < lb:
                    count_s = lb

                ratio = count / count_s
                if ratio < 1:
                    ratio_list_low.append(ratio)
                    q_list_low.append(q)
                    r_list_low.append(r)
                else:
                    ratio_list_high.append(ratio)
                    q_list_high.append(q)
                    r_list_high.append(r)

    df_ratio_low['q'] = q_list_low
    df_ratio_low['r'] = r_list_low
    df_ratio_low['ratio'] = ratio_list_low

    df_ratio_high['q'] = q_list_high
    df_ratio_high['r'] = r_list_high
    df_ratio_high['ratio'] = ratio_list_high

    return df_ratio

def normalize_df(df):
    total_count = df['counts'].sum()
    df['counts'] = df['counts'] / total_count
    return df

def plot_ratio(x, y, names, xs, ys, names_s, lb):

    ratio = (max(y)-min(y)) / (max(x)-min(x))
    size = 20
    # (max(y)-min(y)) / (max(x)-min(x))

    # print(ratio)
    p = figure(title="Hexbin for 500 points", match_aspect=True, aspect_scale = ratio, x_range = [min(x), max(x)],y_range=[min(y),max(y)],
               tools="wheel_zoom,reset", background_fill_color= '#D3D3D3')
    p.grid.visible = False

    alpha = 0.8
    # r, bins = p.hexbin(x, y, size=0.5, hover_color="pink", hover_alpha=0.8)

    df1 = hex_to_bin(x, y, size, ratio)
    df2 = hex_to_bin(xs, ys, size, ratio)
    # df_m2 = hex_to_bin(a, b, size=0.2)

    df1['names'] = names
    df2['names'] = names_s
    # df_m2['names'] = names?

    df1 = df1.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    df2 = df2.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')

    df1 = transform_df(df1)
    df2 = transform_df(df2)

    # normalize counts against total counts
    df1_normalized = normalize_df(df1)
    df2_normalized = normalize_df(df2)

    # calculate ratio's
    df_ratio_low, df_ratio_high = calculate_ratios(df1_normalized, df2_normalized, lb)

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_ratio_low,aspect_scale=ratio,
               fill_color=linear_cmap('ratio', 'Blues9', 0,1))

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_ratio_high,aspect_scale=ratio,
               fill_color=linear_cmap('ratio', 'Blues9', 1, max(df_ratio_high.ratio)))

    hover = HoverTool(tooltips=[("ratio", "@ratio")])
    p.add_tools(hover)
    #
    output_file("plots/compare_ratios.html")

    show(p)

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=False, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    parser.add_argument('-i2', required=False, metavar='input_file_2', dest='input_file_2', help='[i] to select second input file from the results folder')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    file_1 = args.input_file
    file_2 = args.input_file_2
    table_1 = import_table(file_1) # , term ?
    table_2 = import_table(file_2)

    lower_bound = 50 # below this number ...
    x, y, names = create_array(table_1)
    xs, ys, names_s = create_array(table_2)
    plot_ratio(x, y, names, xs, ys, names_s, lower_bound)

if __name__ == '__main__':
    main()
