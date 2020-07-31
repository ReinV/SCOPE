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
from bokeh.palettes import Reds9
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar
from bokeh.layouts import column, row

def hex_to_bin(x, y, size, aspect_scale, orientation="pointytop"):
    # pd = import_required('pandas','hexbin requires pandas to be installed')

    q, r = cartesian_to_axial(x, y, size, orientation, aspect_scale=aspect_scale)
    df = pd.DataFrame(dict(r=r, q=q))
    return df

def import_table(file):
    print('test')
    table = pd.read_pickle(file, compression=None)
    print(table)
    # f = open(file, 'r')
    # lines = f.readlines()
    #
    # table = dict()
    # for line in lines:
    #     line_to_list = line.split('\t')
    #     id = line_to_list[0]
    #     count = int(line_to_list[1])
    #     tfidf = int(line_to_list[2])
    #     name = line_to_list[3]
    #     mass = float(line_to_list[4])
    #     logP = float(line_to_list[5])
    #     superterms = line_to_list[6]
    #     superterms = superterms.strip()
    #     superterms = superterms.strip("\"")
    #     superterms = superterms.strip('[]')
    #     superterms = superterms.replace('\'','')
    #     superterms = superterms.split(', ')
    #
    #     table[id] = {"Count": count,"tfidf": tfidf, "Name": name, "Mass": mass, "logP": logP, "Superterms": superterms}

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
    df_merged = pd.merge(df, dfs, how='outer', on=['q', 'r'])
    df_merged.loc[:,['counts_x', 'counts_y']] = df_merged.loc[:,['counts_x', 'counts_y']].fillna(lb)
    df_merged.loc[:,['counts_x']] = [lb if counts < lb else counts for counts in df_merged.counts_x]
    df_merged.loc[:,['counts_y']] = [lb if counts < lb else counts for counts in df_merged.counts_y]
    df_merged.loc[:,['counts_x']] = df_merged.counts_x / df_merged.counts_x.sum()
    df_merged.loc[:,['counts_y']] = df_merged.counts_y / df_merged.counts_y.sum()
    df_merged.loc[:,'ratio'] = df_merged.counts_x / df_merged.counts_y
    df_merged.loc[:,'log_ratio'] = np.log(df_merged.ratio)

    minimum = min(df_merged.log_ratio)
    maximum = max(df_merged.log_ratio)

    df_ratio_low = df_merged[df_merged.log_ratio < 0]
    df_ratio_high = df_merged[df_merged.log_ratio >= 0]

    print(minimum, maximum)


    return df_ratio_low, df_ratio_high, minimum, maximum

def normalize_df(df):
    total_count = df['counts'].sum()
    df['counts'] = df['counts'] / total_count
    return df

def plot_ratio(x, y, names, xs, ys, names_s, lb, query_1, query_2):

    title = "Hexbin plot comparing %s with %s" % (query_1, query_2)

    max_x = max(np.append(x, xs))
    min_x = min(np.append(x, xs))
    max_y = max(np.append(y, ys))
    min_y = min(np.append(y, ys))

    size = 10
    ratio = (max_y - min_y) / (max_x - min_x)

    # print(ratio)
    p = figure(title=title, match_aspect=True, aspect_scale = ratio, x_range = [min_x, max_x],y_range=[min_y, max_y],
               tools="wheel_zoom,reset, save", background_fill_color= '#D3D3D3')
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
    # df1_normalized = normalize_df(df1)
    # df2_normalized = normalize_df(df2)

    # calculate ratio's
    df_ratio_low, df_ratio_high, minimum, maximum = calculate_ratios(df1, df2, lb)
    extreme = max(abs(minimum), maximum)
    Reds9.reverse()

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_ratio_low,aspect_scale=ratio,
               fill_color=linear_cmap('log_ratio', 'Blues9', -extreme, 0))

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_ratio_high,aspect_scale=ratio,
               fill_color=linear_cmap('log_ratio', Reds9, 0, extreme))

    hover = HoverTool(tooltips=[("log_ratio", "@log_ratio")])
    p.add_tools(hover)
    #

    from bokeh.palettes import Blues7 as blue
    from bokeh.palettes import Reds7 as red

    red.reverse()
    pbr = blue + ['#FFFFFF'] + Reds9

    # color bar
    color_mapper = LinearColorMapper(palette=pbr, low=-5, high=5)

    color_bar = ColorBar(color_mapper=color_mapper, ticker=BasicTicker(),major_label_overrides={4: 'More '+str(query_1)+' counts',-4:"More "+str(query_2)+" counts"},major_label_text_align='left',
                         label_standoff=6, border_line_color=None, location=(0,0))

    # color bar is added to a dummy figure, to prevent the shrinking of the original plot
    dummy = figure(
               toolbar_location=None,
               min_border=0,
               outline_line_color=None)
    dummy.add_layout(color_bar, 'left')
    dummy.title.align = 'center'
    dummy.title.text_font_size = '10pt'

    output_file("plots/compare_ratios.html")

    layout = row(p, dummy)

    show(layout)


def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=False, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    parser.add_argument('-i2', required=False, metavar='input_file_2', dest='input_file_2', help='[i] to select second input file from the results folder')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parser()
    print('wtf')
    file_1 = args.input_file
    file_2 = args.input_file_2
    table_1 = import_table(file_1) # , term ?
    table_2 = import_table(file_2)
    #
    # LOWER_BOUND = 10 # below this number ...
    # x, y, names = create_array(table_1)
    # xs, ys, names_s = create_array(table_2)
    #
    # query_1 = file_1.split('\\')[1].split('_')[0]
    # query_2 = file_2.split('\\')[1].split('_')[0]
    #
    # plot_ratio(x, y, names, xs, ys, names_s, LOWER_BOUND, query_1, query_2)

if __name__ == '__main__':
    main()
