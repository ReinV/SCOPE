#!/usr/bin/python

import numpy as np
import pandas as pd
from math import sqrt, log, exp
import argparse

from collections import Counter
from bokeh.io import output_file, show
from bokeh.models import HoverTool
from bokeh.models import LinearColorMapper, BasicTicker, ColorBar
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import hexbin
from bokeh.util.hex import cartesian_to_axial
from bokeh.models import CustomJS
from bokeh.layouts import row

from bokeh.palettes import Blues7 as blue
from bokeh.palettes import Reds7 as red
from datetime import datetime
from collections import Counter

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
        superterms = line_to_list[6]
        superterms = superterms.strip()
        superterms = superterms.strip("\"")
        superterms = superterms.strip('[]')
        superterms = superterms.replace('\'','')
        superterms = superterms.split(', ')

        table[id] = {"Count": count,"tfidf": tfidf, "Name": name, "Mass": mass, "logP": logP, "Superterms": superterms}

    return table

def create_array(table, normalisation=False):
    '''
    This function recieves the table dictionary and returns numpy arrays with the mass, logP and names of the chemicals.
    These values and names are repeated corresponding to the count value.
    '''
    x = []
    y = []
    names = []
    for id in table.keys():
        if normalisation:
            count = taboe[id]["TFIDF"]
        else:
            count = table[id]["Count"]
        mass = table[id]["Mass"]
        logP = table[id]["logP"]
        name = table[id]["Name"]
        x.extend([logP]*count)
        y.extend([mass]* count)
        names.extend([name]*count)
    return np.asarray(x), np.asarray(y), names

def hex_to_bin(x, y, size, aspect_scale, orientation):
    '''
    This function recieves x, y coordinates arrays and turns these into q, r hexagon coordinates returned in a pandas dataframe.
    '''
    q, r = cartesian_to_axial(x, y, size, orientation, aspect_scale=aspect_scale)
    df = pd.DataFrame(dict(r=r, q=q))
    return df

def count_chemicals(list_names):
    '''
    This function recieves a list with the names of chemicals (possibily many dublicates) in a hexagon (row in a dataframe).
    It returns a dictionary with the unique chemical names as keys and their count as values.
    '''
    uniques = list(set(list_names))
    chemical_to_count = dict()
    for chemical in uniques:
        count = 0
        for name in list_names:
            if chemical == name:
                count += 1
        chemical_to_count[chemical] = count

    return chemical_to_count

def transform_df(df):
    '''
    This function recieves a pandas dataframe with q, r coordinates and names as columns.
    It returns a dataframe with the same coordinates but with a column counts (of all chemicals in that hexagon) and a column names with a dictionary of the chemical names and their individual count in that hexagon.
    Columns are added to the dataframe as lists.
    '''
    counts_list = []
    names_list = []
    for row in df['names']:
        # total count of chemicals
        counts_list.append(len(row))
        # add dictionary with chemical to count
        names_list.append(count_chemicals(row))
    df['counts'] = counts_list
    df['names'] = names_list
    return df

def normalise_df_test(df1,df2,gridcount):

    factor = (df1['counts'].sum()) / (df2['counts'].sum())
    df2['ncounts'] = df2['counts'] * factor
    print(df2['ncounts'].sum())
    print(df1['counts'].sum())
    # df2['ncounts'] += (1-factor)


    # # print(df1[1000:1010])
    # df1['counts'] -= 0.5
    # df2['counts'] -= 0.5
    # # print(df1[1000:1010])
    # df1['counts'] = df1['counts'] / df1['counts'].sum()
    # df2['counts'] = df2['counts'] / df2['counts'].sum()
    # # print(df1[1000:1010])
    return df1, df2

def normalise_df(df1, df2, gridcount):
    '''
    This function recieves a dataframe and returns the same dataframe with normalised counts.
    Counts per row are normalised by dividing with the total count of the whole column.
    '''
    total_count_1 = df1['counts'].sum()-gridcount
    total_count_2 = df2['counts'].sum()-gridcount
    factor = total_count_1 / total_count_2
    df2['ncounts'] = df2['counts'] * factor

    df2['ncounts'] += (1-factor)

    return df1, df2

def calculate_difference(names1, names2):
    '''
    This function recieves two dictionaries of chemical names.
    It returns a dictionary with chemical names in those dictionaries and the absolute differences in count between those dictionaries.
    '''
    chemical_to_difference = dict()

    # go through chemicals in first dictionary
    for key in names1:
        count1 = names1[key]
        # check if chemical is in second dictionary. If so, calculate the absolute difference
        try:
            count2 = names2[key]
            difference = abs(count1 - count2)

        # if chemical is not in second dictionary, the difference is the count in the first dictionary
        except:
            difference = count1
        chemical_to_difference[key] = difference

    # go through chemicals in second dictionary
    for key in names2:
        count2 = names2[key]
        # check if chemical is already evaluated
        try:
            chemical_to_difference[key]
        # if not, add it to the dictionary. No need to check for chemical in the first dictionary because these are already evaluated in first for loop
        except:
            difference = count2
            chemical_to_difference[key] = difference

    return chemical_to_difference

def calculate_ratio_test(df1, df2):
    q_list = []
    r_list = []
    ratio_list = []
    for index, row in df1.iterrows():
        q = row['q']
        r = row['r']
        count1 = row['counts']
        count2 = df2[(df2.q==q) & (df2.r==r)]['ncounts'].iloc[0]
        ratio = count1 / count2
        q_list.append(q)
        r_list.append(r)
        ratio_list.append(ratio)

    df_ratio = pd.DataFrame()
    df_ratio['q'] = q_list
    df_ratio['r'] = r_list
    df_ratio['ratio'] = ratio_list

    return df_ratio

def calculate_ratios(df1, df2):
    '''
    This function recieves two dataframes and calculates the (normalised) count ratio between the dataframes for all q, r coordinates.
    A lowerbound value is applied to normalised counts lower than the lower bound (to prevent extremly high ratio's),
    and to q,r coordinates that do exist in one dataframe but not in the other (so to prevent dividing by zero).
    '''
    lb = 1/10000

    # new dataframe
    df_ratio = pd.DataFrame()

    # columns for the new dataframe
    q_list = []
    r_list = []
    ratio_list = []
    differences = []

    # go through all rows (hexagons) in first dataframe
    for index, row in df1.iterrows():
        q = row['q']
        r = row['r']
        count1 = row['counts']
        names1 = row['names']
        if count1 < lb:
            count1 = lb # apply lower bound ..

        # check if row exists in second dataframe
        try:
            count2 = df2[(df2.q==q) & (df2.r==r)]['ncounts'].iloc[0]
            names2 = df2[(df2.q==q) & (df2.r==r)]['names'].iloc[0]
            if count2 < lb:
                count2 = lb # apply lower bound ..

        # if not, apply lowerbound
        except:
            count2 = lb # apply lower bound ..
            names2 = dict()

        # add ratio to lists that will be added to a new dataframe as columns
        ratio = count1/count2
        ratio_list.append(ratio)
        q_list.append(q)
        r_list.append(r)

        # calculate the differences between the counts of chemicals of the two dataframes in a hexagon
        chemical_to_difference = calculate_difference(names1, names2)
        differences.append(chemical_to_difference)

    # go through the rows of the second dataframe
    for index, row in df2.iterrows():
        q = row['q']
        r = row['r']
        counts2 = row['counts']
        if counts2 < lb:
            counts2 = lb # apply lowe bound

        #check if coordinates exists in df1, if so, this row can be skipped
        try:
            counts1 = df1[(df1.q==q) & (df1.r==r)]['counts'].iloc[0]

        # if not, apply lower bound and add values to the list
        except:
            count1 = lb # apply lower bound
            names1 = dict()
            ratio = count1/count2
            ratio_list.append(ratio)
            q_list.append(q)
            r_list.append(r)

            chemical_to_difference = calculate_difference(names1, names2)
            differences.append(chemical_to_difference)

    # add lists as columns to the new dataframe
    df_ratio['q'] = q_list
    df_ratio['r'] = r_list
    df_ratio['ratio'] = ratio_list
    # df_ratio['dif'] = differences

    return df_ratio

def find_discriminating_chemicals(df):
    '''
    This function recieves a dataframe, and adds to this dataframe 5 columns.
    These columns consist of the xth most discriminating chemical, determined by the absolute differences already calculated and saved in the dictionaries in the 'dif' column.
    '''
    # make columns
    columns = [[] for _ in range(5)]

    for row in df['dif']:
        # dictionary with values
        chemical_to_difference = row
        # function that returns a list of tuples (key, value)
        k = Counter(chemical_to_difference)
        # function that returns the 5 tuples based on the highest values
        high = k.most_common(5)
        # get the chemicals + counts and put them in the column lists
        for i in range(5):
            if i < len(high):
                name = high[i]
            # sometimes there are less than 5 chemicals in a hexagon, the add this tuple to the column
            else:
                name = ('-','-')
            columns[i].append(name)

    #add lists to df as columns
    for i in range(5):
        column_name = "names"+str(i+1)
        df[column_name] = columns[i]

    return df

def get_blur(x,y):
    '''
    This function recieves x, y values and sigma x, sigma y values and returns the calculated blur value.
    See https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    '''
    sigma_x = 2
    sigma_y = 1
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def add_gaussian_blur(df):
    '''
    This function adds gaussian blur to the plot by checking every possible coordinate with a gaussian kernel and calculating a new value.
    The kernel has a center, inside ring and outside ring with attached weights.
    '''

    # construct kernel
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

    # add blur column
    df['blur'] = df['log']

    # make normal python dict
    coords_to_count = dict()
    coords_to_blur = dict()

    for index, row in df.iterrows():
        q = row['q']
        r = row['r']
        log_ratio = row['log']
        coords_to_blur[(q,r)] = {}

        # add dataframe to a python dictionary
        coords_to_count[(q,r)] = log_ratio
        for i in range(len(df.columns)):
            column_name = df.columns[i]
            column_value = row[i]
            coords_to_blur[(q,r)][column_name] = column_value

    # go through all coordinates
    for coord in coords_to_count.keys():
        q = coord[0]
        r = coord[1]
        log_ratio = coords_to_count[coord]
        # for every coordinate, apply kernel (to surrounding hexagons)
        for key in kernel.keys():
            q_new = q + key[0]
            r_new = r + key[1]
            factor = kernel[key]
            blur_new = log_ratio*factor
            # check if surrounding hexagon exists and apply blur
            try:
                coords_to_blur[(q_new, r_new)]['blur'] += blur_new
            except:
                coords_to_blur[(q_new, r_new)] = {'q': q_new, 'r': r_new, 'ratio': 1, 'log': 0, 'blur': blur_new}

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

def make_grid(x_min, x_max, y_min, y_max, w, h):
    coords_x_1 = np.arange(x_min-1, x_max+1, 1.5*w)
    coords_x_2 = np.arange(x_min+w-1, x_max+w+1, 1.5*w)
    coords_y = np.arange(y_min, y_max, h)
    x_list = []
    y_list = []
    for element in coords_x_1:
        x_list = np.append(x_list, np.repeat(element, len(coords_y)) )
        y_list = np.append(y_list, coords_y)
    for element in coords_x_2:
        x_list = np.append(x_list, np.repeat(element, len(coords_y)))
        y_list = np.append(y_list, coords_y)
    return x_list, y_list

def plot_ratio(x1, y1, names1, x2, y2, names2, term1, term2, lb):
    '''
    This function recieves the (repeated) x,y coördinates, names of the chemicals, and the terms of two different searches.
    Counts are determined for the hexagons. After normalisation for total counts, a log(ratio) is calculated for every hexagon.
    Then, a blur is applied to the log(ratio), but the original log(ratio) is shown in the tooltip.
    The hexagon figure is saved in a html file.
    '''

    length = len(x1) + len(x2)
    title = 'Hexbin plot for ' + str(length) + ' annotated chemicals with ratios of counts from ' + str(term1) + '/' + str(term2)

    x_min = min(min(x1),min(x2))
    x_max = max(max(x1),max(x2))
    y_min = min(min(y1),min(y2))
    y_max = max(max(y1),max(y2))

    # hexagon properties
    orientation = 'flattop'
    size = 20
    ratio = (y_max-y_min) / (x_max-x_min )
    if orientation == 'flattop':
        size = size / ratio
    h = (sqrt(3) * size) * ratio
    w = (2 * size)

    # boundaries plot / grid

    x_grid, y_grid = make_grid(x_min, x_max, y_min, y_max, w, h)

    # make figure
    p = figure(title=title, x_range = [x_min-0.5, x_max+0.5],y_range=[0-(h/2),y_max+100-h],
    tools="wheel_zoom,reset", background_fill_color= '#FFFFFF')
    #'#D3D3D3'
    p.grid.visible = False

    x1 = np.append(x1, x_grid)
    y1 = np.append(y1, y_grid)
    x2 = np.append(x2, x_grid)
    y2 = np.append(y2, y_grid)

    df1 = hexbin(x1, y1, size, orientation, ratio)
    df2 = hexbin(x2, y2, size, orientation, ratio)

    df1, df2 = normalise_df_test(df1, df2, len(x_grid))

    df3 = calculate_ratio_test(df1, df2)
    df3['log'] = np.log(df3['ratio'])

    print(df3['log'].sum())
    # df3['blur'] = df3['log']
    df3 = add_gaussian_blur(df3)

    print(df3['log'].sum())

    highest_abs_value = max(max(df3.blur),abs(min(df3.blur)))

    df_low = df3[df3['blur'] < 0]
    df_high = df3[df3['blur'] > 0]
    df_eq = df3[df3['blur'] == 0]

    red.reverse()
    pbr = blue + ['#FFFFFF'] + red


    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_low, aspect_scale=ratio, orientation=orientation, hover_color="#39FF14",
               fill_color=linear_cmap('blur', blue, -highest_abs_value, 0))

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_high,aspect_scale=ratio, orientation=orientation, hover_color="#39FF14",
               fill_color=linear_cmap('blur', red, 0, highest_abs_value))

    p.hex_tile(q="q", r="r", size=size, line_color=None, source=df_eq,aspect_scale=ratio, orientation=orientation, hover_color="#39FF14",
               fill_color='#FFFFFF')

    #hover tool
    hover = HoverTool(tooltips=[('blur', '@blur')])
    # hover = HoverTool(tooltips=[("log(ratio)", "@log"),("#1", "@names1"),("#2", "@names2"),("#3", "@names3"),("#4", "@names4"),("#5", "@names5")],callback=callback_hover, show_arrow=False)
    p.add_tools(hover)

    #output
    file_name = "plots/compare_ratios_" + str(term1) + "_" + str(term2) + ".html"
    output_file(file_name)

    # r = row(p,dummy)

    show(p)

    # p.hex_tile(q="q", r="r", size=size, line_color=None, source=df3, aspect_scale=ratio, orientation=orientation, hover_color='#39FF14',
    # fill_color=linear_cmap('ratio', 'Viridis256', 0, max(df3.ratio) ))
    #
    # print(len(df_low))
    # print(len(df_high))
    # print(len(df_eq))


    # print(df_low['log'].sum())
    # print(df_high['log'].sum())
    # get coordinates for every x,y point
    # df1 = hexbin(x3, y3, size, orientation, ratio)
    # df1 = hex_to_bin(x1, y1, size, ratio, orientation)
    # df2 = hex_to_bin(x2, y2, size, ratio, orientation)
    #
    # # add the names to the dataframe
    # df1['names'] = names1.append(['']*len(x3))
    # df2['names'] = names2.append(['']*len(x3))
    #
    # # group by coordinates for hexagons
    # df1 = df1.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    # df2 = df2.groupby(['q', 'r'])['names'].apply(list).reset_index(name='names')
    #
    # # count chemicals per hexagon, also add column with a list of counts of individual chemicals # (SHOULD THESE BE NORMALISED?)
    # df1 = transform_df(df1)
    # df2 = transform_df(df2)
    # #
    # # print(len(df1))
    # # print(len(df2))
    # #
    # # normalize counts against total counts
    # df1_normalized = normalize_df(df1)
    # df2_normalized = normalize_df(df2)
    #
    # # print(df1_normalized['counts'].sum())
    # # print(df2_normalized['counts'].sum())
    # #
    # # # calculate ratio's and find most discriminating chemicals (+ difference? in counts)
    # df_ratio = calculate_ratio_test(df1_normalized, df2_normalized)
    #
    # df_ratio_low = df_ratio[df_ratio['ratio'] < 1]
    # df_ratio_eq = df_ratio[df_ratio['ratio'] == 1]
    # df_ratio_high = df_ratio[df_ratio['ratio'] > 1]
    #
    # df_ratio_low['log'] = np.log(df_ratio_low['ratio'])
    # df_ratio_high['log'] = np.log(df_ratio_high['ratio'])

    # print(len(df_ratio_low))
    # print(len(df_ratio_high))
    # print(df_ratio_high['log'].sum())
    # print(df_ratio_low['log'].sum())
    # print(df_ratio[0:10])
    #
    # # log transformation
    # df_ratio['log'] = np.log(df_ratio['ratio'])
    # df_ratio['blur'] = np.log(df_ratio['ratio'])
    # # finding most discriminating chemicals per hexagon
    # # df_ratio = find_discriminating_chemicals(df_ratio)
    # # df_ratio = df_ratio.drop(columns="dif")
    #
    # print(df_ratio.columns)
    # # add blur to log ratio's
    # df_ratio = add_gaussian_blur(df_ratio)
    # print(df_ratio.columns)
    #
    # # make two dataframes to create two glyphs
    # df_ratio_low = df_ratio[df_ratio['log'] < 0]
    # df_ratio_zero = df_ratio[df_ratio['log'] == 0]
    # df_ratio_high = df_ratio[df_ratio['log'] > 0]
    #
    # print(len(df_ratio_low))
    # print(len(df_ratio_high ))
    # print(len(df_ratio_zero))
    # # blue.append('#FFFFFF')
    # # red.append('#FFFFFF')

    # # color bar
    # color_mapper = LinearColorMapper(palette=pbr, low=-5, high=5)
    #
    # color_bar = ColorBar(color_mapper=color_mapper, ticker=BasicTicker(),major_label_overrides={4: 'More '+str(term1)+' counts',-4:"More "+str(term2)+" counts"},major_label_text_align='left',
    #                      label_standoff=6, border_line_color=None, location=(0,0))
    #
    # # color bar is added to a dummy figure, to prevent the shrinking of the original plot
    # dummy = figure(
    #            toolbar_location=None,
    #            min_border=0,
    #            outline_line_color=None)
    # dummy.add_layout(color_bar, 'left')
    # dummy.title.align = 'center'
    # dummy.title.text_font_size = '10pt'


    #hover position
    callback_hover = CustomJS(code="""
        var tooltips = document.getElementsByClassName("bk-tooltip");
        for (var i = 0, len = tooltips.length; i < len; i ++) {
            tooltips[i].style.top = ""; // unset what bokeh.js sets
            tooltips[i].style.left = "";
            tooltips[i].style.bottom = "280px";
            tooltips[i].style.left = "700px";
            tooltips[i].style.width = "300px";
        }
        """)

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=False, metavar='input_file', dest='input_file', help='[i] to select input file from the results folder')
    parser.add_argument('-i2', required=False, metavar='input_file_2', dest='input_file_2', help='[i] to select second input file from the results folder')
    arguments = parser.parse_args()
    return arguments

def main():
    startTime = datetime.now()
    args = parser()
    file_1 = args.input_file
    file_2 = args.input_file_2

    term_1 = file_1.split('\\')[1].split('_')[0]
    term_2 = file_2.split('\\')[1].split('_')[0]

    table_1 = import_table(file_1) # , term ?
    table_2 = import_table(file_2)

    lower_bound = 50 # below this number ...
    x1, y1, names1 = create_array(table_1)
    x2, y2, names2 = create_array(table_2)
    plot_ratio(x1, y1, names1, x2, y2, names2, term_1, term_2, lower_bound)

    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
