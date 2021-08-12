#!/usr/bin/python

import numpy as np
import pandas as pd
import math
import argparse
import re

from collections import Counter
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import hexbin
from bokeh.util.hex import cartesian_to_axial
from bokeh.palettes import Blues9
from bokeh.palettes import Reds9
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar
from bokeh.layouts import column, row

def create_data_source(table, size, orientation, ratio):
    '''
    This function recieves x and y vectors and returns a dataframe with hexiconal coordinates + counts
    '''
    # Create array with mass and logP values
    x, y = create_array(table)
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)
    # create dataframe with hexagonal coordinates as rows (row = hexagon)
    df = table.reset_index()
    df.loc[:,"q"] = q
    df.loc[:,"r"] = r
    df = df.drop(['ChEBI', 'Class', 'Names'], axis=1)
    df = df.groupby(['q', 'r']).agg({'Count': 'sum', 'TFIDF': 'sum'}).reset_index()
    return df

def import_table(file):
    '''
    This function imports a pickle file (from the 'table' folder) and returns a dataframe
    '''
    table = pd.read_pickle(file)
    return table

def create_array(table):
    '''
    This function recieves a dataframe with logP and Mass values for every ChEBI identifier.
    It returns two numpy arrays: for mass and logP.
    '''
    # create lists
    x = [float(logP) for logP in table.logP]
    y = [float(mass) for mass in table.Mass]

    return np.asarray(x), np.asarray(y)

def get_ratio(x, y, total_x, total_y, lb, normalize):
    '''
    This function is used to calculate the ratio between counts of two dataframes in one hexagon.
    If the normalize argument is given, then both counts in the hexaon will be normalized for the total amount of counts found in the dataframes.
    '''
    if normalize:
        factor = total_x / total_y
        if total_x < total_y:
            if x > lb:
                x = x / factor
        else:
            if y > lb:
                y = y * factor
    ratio = x / y
    return ratio

def calculate_ratios(df1, df2, lb, normalize):
    '''
    This function recieves two dataframes, merges these together and calculates the count ratio's for every hexagon.
    The ratio per hexagon will be represented in the plot with red (high ratio), blue (low ratio), and white (equal) colors.
    Two dataframes for low and high ratio levels are returned.
    '''

    # Merge dataframes
    df_merged = pd.merge(df1, df2, how='outer', on=['q', 'r'])

    # Apply lower bound to every hexagon
    df_merged.loc[:,['Count_x', 'Count_y']] = df_merged.loc[:,['Count_x', 'Count_y']].fillna(lb)
    df_merged.loc[df_merged.Count_x < lb, "Count_x"] = lb
    df_merged.loc[df_merged.Count_y < lb, "Count_y"] = lb

    # Get total counts
    total_count_x = df_merged.Count_x.sum()
    total_count_y = df_merged.Count_y.sum()

    # Calculate ratio's per hexagon
    df_merged.loc[:,'ratio'] = [get_ratio(x, y, total_count_x, total_count_y, lb, normalize) for x, y in zip(df_merged.Count_x, df_merged.Count_y)]
    df_merged.loc[:,'log_ratio'] = np.log(df_merged.ratio)

    # Divide dataframe into low and high ratio dataframe (for plot colors)
    minimum = min(df_merged.log_ratio)
    maximum = max(df_merged.log_ratio)
    df_ratio_low = df_merged[df_merged.log_ratio < 0]
    df_ratio_high = df_merged[df_merged.log_ratio >= 0]

    return df_ratio_low, df_ratio_high, minimum, maximum

def normalize_df(df):
    '''
    This function applies normalization tot the total counts in a dataframe
    '''
    total_count = df['counts'].sum()
    df['counts'] = df['counts'] / total_count
    return df

def find_max_values(table_1, table_2):
    '''
    This function recieves the input tables, creates x and y arrays and finds the min and max values for the plot.
    '''
    # Find max values to calculate ratio for the plot
    x, y = create_array(table_1)
    xs, ys = create_array(table_2)
    max_x = max(np.append(x, xs))
    min_x = min(np.append(x, xs))
    max_y = max(np.append(y, ys))
    min_y = min(np.append(y, ys))
    ratio = (max_y - min_y) / (max_x - min_x)
    return max_x, min_x, max_y, min_y, ratio

def plot_ratio(table_1, table_2, query_1, query_2, normalize):
    # Define constants
    SIZE = 10
    LOWER_BOUND = 2
    ORIENTATION = "pointytop"

    # Plot variables
    title = "Hexbin plot comparing %s with %s" % (query_1, query_2)
    max_x, min_x, max_y, min_y, ratio = find_max_values(table_1, table_2)

    # Create figure
    p = figure(title=title, match_aspect=True, aspect_scale = ratio, x_range = [min_x, max_x],y_range=[min_y, max_y],
               tools="wheel_zoom,reset, save", background_fill_color= '#D3D3D3')
    p.grid.visible = False

    # Create hexagonal dataframes (turn LogP and Mass values into hexagonal coordinates)
    df1 = create_data_source(table_1, SIZE, ORIENTATION, ratio)
    df2 = create_data_source(table_2, SIZE, ORIENTATION, ratio)

    # Calculate ratio's for plot dataframes, and min max values for
    df_ratio_low, df_ratio_high, minimum, maximum = calculate_ratios(df1, df2, LOWER_BOUND, normalize)
    extreme = max(abs(minimum), maximum)

    # Reverse list of colors for color bar
    red_reversed = list(Reds9)
    red_reversed.reverse()

    # Create blue and red hex tiles for low and high ratio's.
    p.hex_tile(q="q", r="r", size=SIZE, line_color=None, source=df_ratio_low,aspect_scale=ratio,
               fill_color=linear_cmap('log_ratio', 'Blues9', -extreme, 0))

    p.hex_tile(q="q", r="r", size=SIZE, line_color=None, source=df_ratio_high,aspect_scale=ratio,
               fill_color=linear_cmap('log_ratio', red_reversed, 0, extreme))

    hover = HoverTool(tooltips=[("log_ratio", "@log_ratio")])
    p.add_tools(hover)

    pbr = list(Blues9) + ['#FFFFFF'] + red_reversed

    # Color bar
    color_mapper = LinearColorMapper(palette=pbr, low=-5, high=5)
    color_bar = ColorBar(color_mapper=color_mapper, ticker=BasicTicker(),major_label_overrides={4: 'More '+str(query_1)+' counts',-4:"More "+str(query_2)+" counts"},major_label_text_align='left',
                         label_standoff=6, border_line_color=None, location=(0,0))

    # Color bar is added to a dummy figure, to prevent the shrinking of the original plot
    dummy = figure(
               toolbar_location=None,
               min_border=0,
               outline_line_color=None)
    dummy.add_layout(color_bar, 'left')
    dummy.title.align = 'center'
    dummy.title.text_font_size = '10pt'

    # Output
    output_file("plots/%s_vs_%s.html" % (query_1, query_2))
    layout = row(p, dummy)
    show(layout)


def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=False, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    parser.add_argument('-i2', required=False, metavar='input_file_2', dest='input_file_2', help='[i] to select second input file from the results folder')
    parser.add_argument('-n', default=False, action='store_true', help='[n] to select normalization')
    arguments = parser.parse_args()
    return arguments

def main():
    # Arguments
    args = parser()
    file_1 = args.input_file
    file_2 = args.input_file_2
    normalize = args.n

    # Import input files
    table_1 = import_table(file_1) # , term ?
    table_2 = import_table(file_2)

    # Get query names from input (Windows/Linux)
    query_1 = re.split(r'[\\/]',file_1)[1].split('_')[0]
    query_2 = re.split(r'[\\/]',file_2)[1].split('_')[0]

    # Plot!
    plot_ratio(table_1, table_2, query_1, query_2, normalize)

if __name__ == '__main__':
    main()
