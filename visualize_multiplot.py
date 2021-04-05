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

# BOKEH
from bokeh import events
from bokeh.io import output_file, show
from bokeh.models import CustomJS, HoverTool, ColumnDataSource, Slider, CheckboxGroup, RadioGroup, Button, MultiSelect
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.transform import log_cmap
from bokeh.util.hex import axial_to_cartesian
from bokeh.util.hex import cartesian_to_axial
from bokeh.layouts import column, row
from bokeh.palettes import Viridis256, Greys256

def import_table(file):
    '''
    This function imports the pkl files from the tables forlder, and returns a pandas dataframe.
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

def hexbin(df, x, y, size, aspect_scale, orientation):
    '''
    This function recieves x and y coordinate arrays and converts these into q and r hexagon coordinates by calling Bokeh's "cartesian_to_axial" function.
    The q and r coordinates are added to the dataframe, and this dataframe is returend.
    '''
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=aspect_scale)

    df.loc[:'q'] = q
    df.loc[:'r'] = r

    return df

def add_tooltip_columns(df, table):
    '''
    For every hexagon, a tooltip will be created that will be shown when the user hovers with the mouse over the hexagon.
    The tooltip will show the 3 most frequent ChEBI identifiers and additional information.

    This function recieves:
    - a "df" dataframe (with hexagonal coordinates) that will be the source for the multiplot
    - the "table" dataframe with information that is needed for the tooltip such as name, mass, logP for every ChEBI identifier

    In this function, the tooltip information will be added to the original dataframe in additional columns.
    These columns will be used by JavaScript code to display in the tooltip.

    '''
    table = table.drop(['Class', 'logP', 'Mass'], axis=1)
    table = table.reset_index()

    # Define tooltip size
    TOOLTIP_COUNT = 3
    columns = table.columns
    tooltip_columns = {column+str(i):[] for i in range(1, TOOLTIP_COUNT+1) for column in columns}

    # Extract ChEBI identifiers from dataframe
    chebi_ids = [ids if isinstance(ids, list) else [] for ids in df.ChEBI]
    list_for_df = []

    # Use chebi identifiers to look up information in the "table" dataframe
    for ids in chebi_ids:
        rows = table[table.ChEBI.isin(ids)]

        # Sort and select most frequent ChEBI identifiers
        rows = rows.sort_values(by='Count', ascending=False)
        values_nested = rows[0:TOOLTIP_COUNT-1].values.tolist()

        # Unnest information from the table,
        values_unnested = [item for sublist in values_nested for item in sublist]
        while len(values_unnested) < (TOOLTIP_COUNT*len(columns)):
            values_unnested.append("-")

        list_for_df.append(values_unnested)

    df_tooltip = pd.DataFrame(list_for_df, columns = tooltip_columns)

    df = df.join(df_tooltip, how='left')

    return df

def get_blur(x,y,sigma_x,sigma_y):
    '''
    This function recieves x, y values and sigma x, sigma y values and returns the calculated blur value.
    See https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    '''
    return exp(-0.5*(x*x/sigma_x/sigma_x + y*y/sigma_y/sigma_y))

def get_rows(q, r, counts, tfidf, kernel, blur_max, step_size):
    '''
    For every hexagon (=row in the dataframe), counts will be distributed to surrounding hexagons in a gaussian manner (blur).
    This function recieves information for one hexagon so that it can distribute its counts to other hexagons.
    To distribute the counts to other hexagons/rows, we use the "kernel".
    We create new rows for other hexagons, which means we might create many rows with the same hexagonal coorindates.
    Newly created rows will be returned
    '''
    # initiate original row
    rows = [[q, r, counts, tfidf] + [counts for i in np.arange(0, blur_max+step_size, step_size)] + [tfidf for i in np.arange(0, blur_max+step_size, step_size)]]

    # use kernel to distribute counts to other rows
    for coords_new, blur_factors in kernel.items():
        q_new = q + coords_new[0]
        r_new = r + coords_new[1]

        # create new row, use kernels coordinates and calculated blur values
        new_row = [q_new, r_new, 0, 0] + list(map(lambda n: n * counts, blur_factors)) + list(map(lambda n: n*tfidf, blur_factors))
        rows.append(new_row)
    return rows

def construct_kernel(blur_max, step_size):
    '''
    This function recieves the maximum blur and step size to construct the kernel.
    The kernel is a dictionary that uses the coordinates as keys, and the blur values as values.
    The blur values depend on sd_x, so the kernel will return a list of blur values from 0 to blur_max.
    Blur values are calculated by "get_blur()"".
    The kernel will be used by "get_rows()" to distribute counts to surrounding hexagons.
    '''
    coordinates_to_distance = {(-5,2):(7.5, sqrt(3)/2),(-5,3):(7.5, sqrt(3)/2),
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

    kernel = {}
    for key, distance in coordinates_to_distance.items():
        kernel[key] = [0]
        for sd_x in np.arange(step_size, blur_max+step_size, step_size):
            kernel[key].append(get_blur(distance[0], distance[1], sd_x, sd_x/2))

    return kernel

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

    kernel = construct_kernel(blur_max, step_size)

    columns = ['q', 'r', 'Count', 'TFIDF'] + [str(sd_x) for sd_x in np.arange(0, (blur_max+step_size), step_size)] + ['%s_tfidf' % str(sd_x) for sd_x in np.arange(0, (blur_max+step_size), step_size)]

    df_blur = pd.concat([pd.DataFrame(get_rows(q, r, counts, tfidf, kernel, blur_max, step_size),
        columns=columns) for q, r, counts, tfidf in zip(df.q, df.r, df.Count, df.TFIDF)], ignore_index=True)

    df_blur = df_blur.groupby(['q', 'r'], as_index=False).agg(sum)

    df_joined = df_blur.merge(df.loc[:,['q', 'r', 'ChEBI']], on=['q', 'r'], how='outer')

    # aggregation_dictionary = {}
    # for column in df_blur.columns:
    #     if column != 'q' and column != 'r':
    #         aggregation_dictionary[column] = 'sum'
    #
    # df_blur = df_blur.groupby(['q', 'r'], as_index=False).agg(aggregation_dictionary)

    return df_joined

def check_for_ids(id_list, class_id):
    '''
    This functions checks if the Class identifier is in the "id_list".
    This ChEBI ID-specific list contains other ChEBI identifers of chemicals that are of higher level hierarchical class.
    '''
    # check = any(class_id == id for id in id_list)
    check = any(str(class_id) == str(id) for id in id_list)
    return check

def create_class_source(table, size, ratio, orientation, class_id):
    '''
    This function finds all chemicals belonging to class "class_id" as defined by the ChEBI ontology.
    It returns a dataframe filtered for these chemicals.
    '''
    if class_id == None:
        # if no class id is given, then class source should be empty
        df = pd.DataFrame().reindex_like(table)

    else:
        table_class_only = table[[check_for_ids(id_list, class_id) for id_list in table["Class"]]]

        # create array
        x, y = create_array(table_class_only)
        q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)
        df = table_class_only.reset_index()
        df.loc[:,"q"] = q
        df.loc[:,"r"] = r

        df = df.groupby(['q', 'r']).agg({'Count': 'sum', 'TFIDF': 'sum', 'ChEBI': list}).reset_index()
        df = df.drop(columns='ChEBI')
    # Create plot source
    source = ColumnDataSource(df)

    return source

def create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE):
    '''
    This function recieves a query-specific "table" dataframe, and constructs a source for the hexagonal plot.
    The plot source will be created depending on plot-specific values (size, ratio, orientation).
    Additionally, a gaussian blur will be applied for different values of sd(x) (BLUR_MAX, BLUR_STEP_SIZE)
    '''
    # create array with mass and logP values
    x, y = create_array(table)

    # create hexagonal coordinates from mass and logP values
    q, r = cartesian_to_axial(x, y, size, orientation=orientation, aspect_scale=ratio)

    # create dataframe with hexagonal coordinates as rows (row = hexagon)
    df = table.reset_index()
    df.loc[:,"q"] = q
    df.loc[:,"r"] = r

    # sum rows with identical coordinates together, add blur, and add tooltip information
    df = df.groupby(['q', 'r']).agg({'Count': 'sum', 'TFIDF': 'sum', 'ChEBI': list}).reset_index()
    df = add_gaussian_blur(df, BLUR_MAX, BLUR_STEP_SIZE)
    df = add_tooltip_columns(df, table)
    df = df.drop(columns='ChEBI')

    # plot title and source
    title = 'Hexbin plot for '+str(len(x))+' annotated chemicals with query '+str(term)
    source = ColumnDataSource(df)
    return source, title


def return_JS_code(widget):
    '''
    This function recieves a string that indicates the widget.
    It returns a JavaScript callback code corresponding to the widget.
    '''
    if widget == 'multi_select':
        code = """
            var source_data = source.data;
            var source_class_data = source_class.data;
            var term_to_source = term_to_source;
            var term = cb_obj.value[0];
            var f = slider1.value
            var sd_x = slider2.value;
            var p = p;
            var mapper = mapper;
            var class_hex = class_hex;
            var checkbox_class = checkbox_class
            var term_to_class = term_to_class;


            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf

            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // select new_data
            var new_data = term_to_source[term]['source'].data
            var new_class = term_to_class[term]['source'].data

            // replace current data
            for (var key in new_data) {
                source_data[key] = [];
                for (i=0;i<new_data[key].length;i++) {
                    source_data[key].push(new_data[key][i]);
                }
            }

            // class

            class_hex.visible = false

            for (var key in new_class){
                source_class_data[key] = [];
                for (i=0;i<new_class[key].length;i++) {
                    source_class_data[key].push(new_class[key][i])
                }
            }

            class_hex.source = source_class_data
            checkbox_class.active = []

            // new title
            var title = term_to_source[term]['title']
            p.title.text = title

            // apply blur and saturation

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source_class.change.emit();
            source.change.emit();
        """
    elif widget == 'slider2':
        code = """
            var source_data = source.data;
            var f = slider1.value
            var mapper = mapper;
            var checkbox = checkbox;
            var sd_x = cb_obj.value;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf

            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply blur and saturation

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill_color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
        """
    elif widget == 'tooltips':
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
              <col width="80">
              <col width="305">
              <tr>
                <th>Total Counts</th>
                <th>@Count</th>
                <th>@TFIDF</th>
                <th>($x, $y)</th>
              </tr>
              <tr style="color: #fff; background: black;">
                <th>Chebi ID</th>
                <th>Count</th>
                <th>TFIDF</th>
                <th>Name</th>
              </tr>
              <tr>
                <th>@ChEBI1</th>
                <th>@Count1</th>
                <th>@TFIDF1</th>
                <th>@Names1</th>
              </tr>
              <tr>
                <th>@ChEBI2</th>
                <th>@Count2</th>
                <th>@TFIDF2</th>
                <th>@Names2</th>
              </tr>
              <tr>
                <th>@ChEBI3</th>
                <th>@Count3</th>
                <th>@TFIDF3</th>
                <th>@Names3</th>
              </tr>
            </table>

            """

    elif widget == 'slider1':
        code = """
            var mapper = mapper
            var source_data = source.data;
            var f = cb_obj.value;
            var checkbox = checkbox;
            var sd_x = slider2.value;
            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf

            if (checkbox.active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply scaling

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }


            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
            source.change.emit();
            """

    elif widget == 'rbg':
        code = """
            var active = cb_obj.active;
            var p = p;
            var Viridis256 = Viridis256;
            var Greys256 = Greys256;
            // var class_hex = class_hex;
            // var term_to_class = term_to_class;
            var term = multi_select.value[0]

            if (active == 1){
            mapper.transform.palette = Greys256
            p.background_fill_color = '#000000'
            // if (term_to_class[term]['show_class']){class_hex.glyph.fill_color = 'red'}
            }

            if (active == 0){
            mapper.transform.palette = Viridis256
            p.background_fill_color = '#440154'
            // if(term_to_class[term]['show_class']){class_hex.glyph.fill_color = 'pink'}
            }

            """

    elif widget == 'checkbox':
        code = """
            var source_data = source.data;
            var active = cb_obj.active
            var f = slider1.value
            var mapper = mapper
            var sd_x = slider2.value;

            if (sd_x % 1 == 0){var sd_x = sd_x.toFixed(1)}
            var sd_x = String(sd_x);

            // check for tfidf

            if (active.length == 1) {
            sd_x = sd_x.concat('_tfidf')
            }

            // apply scaling

            for (var i = 0; i < source_data[sd_x].length; i++) {
                source_data['Count'][i] = Math.pow(source_data[sd_x][i], 1/f)
                }

            // maximum value for fill color
            mapper.transform.high = Math.max.apply(Math, source_data['Count'])

            // apply changes
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
    elif widget == 'stats':
        code = """
            var term_to_stats = term_to_stats
            var term = multi_select.value[0]
            var stats_description = term_to_stats[term]
            var wnd = window.open("about:blank", "", "_blank");
            wnd.document.write(stats_description)
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

def create_stats_description(table):
    '''
    This function recieves the "table" variable and retrieves summary statistcs (mean, std).
    Statistics are returned in html style.
    '''
    total_count = table.Count.sum()
    mean_logP = np.repeat(table.logP, table.Count).values.astype(float).mean()
    mean_mass = np.repeat(table.Mass, table.Count).values.astype(float).mean()
    stdev_logP = np.repeat(table.logP, table.Count).values.astype(float).std()
    stdev_mass = np.repeat(table.Mass, table.Count).values.astype(float).std()
    title = 'Statistics of mass and logP values'
    stats_listed = [title, 'Total amount of chemicals: %d' % total_count , 'LogP mean: %f' % mean_logP, 'LogP standard deviations: %f' % stdev_logP,
    'Mass mean: %f' % mean_mass, 'Mass standard deviation: %f' % stdev_mass]
    stats_description = return_html(stats_listed)
    return stats_description

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

def plot(tables, output_filename, xmin, xmax, ymin, ymax, class_id):
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

    # source for widgets
    term_to_source = dict()
    term_to_class = dict()
    term_to_metadata = dict()
    term_to_stats = dict()

    options = []
    # Loop for plot sources
    for term in tables.keys():
        print('sourcing %s' % term)

        # add term to the options for the multiplot
        options.append((term, term))

        # get table
        table = tables[term]['table']

        # plot sources
        source, title = create_data_source(table, term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
        source_class = create_class_source(table, size, ratio, orientation, class_id)
        stats_description = create_stats_description(table)
        term_to_source[term] = {'source': source, 'title': title}
        term_to_class[term] =  {'source': source_class, 'show_class': True}
        term_to_stats[term] = stats_description

        # metadata
        metadata = return_html(tables[term]['metadata'])
        term_to_metadata[term] = metadata

    # make default souce for plot, this is the first source shown in the plot, and also works like a container. Old data is thrown out and new data is thrown in.
    default_term = list(tables.keys())[0] # pick the first one
    table = tables[default_term]['table']
    source, title = create_data_source(table, default_term, size, ratio, orientation, BLUR_MAX, BLUR_STEP_SIZE)
    p.title.text = title
    metadata = tables[default_term]['metadata']
    metadata = return_html(metadata)

    # color mapper
    mapper = linear_cmap('Count', 'Viridis256', 0, max(source.data['Count']))
    hex = p.hex_tile(q="q", r="r", size=size, line_color=None, source=source, aspect_scale=ratio, orientation=orientation,
           fill_color=mapper)

    if class_id:
        source_class = create_class_source(table, size, ratio, orientation, class_id)
        class_hex = p.hex_tile(q='q', r="r", size=size, line_color=None, source=source_class, aspect_scale=ratio,orientation=orientation,
            fill_color='#ff007f')
        class_hex.visible = False

    # HOVER
    TOOLTIPS = return_JS_code('tooltips')
    code_callback_hover = return_JS_code('hover')
    callback_hover = CustomJS(code=code_callback_hover)
    hover = HoverTool(renderers=[hex], tooltips=TOOLTIPS, callback=callback_hover, show_arrow=False)
    p.add_tools(hover)

    # Widgets
    slider1 = Slider(start=1, end=SATURATION_MAX, value=1, step=SATURATION_STEP_SIZE, title="Saturation", width=100)
    slider2 = Slider(start=0, end=BLUR_MAX, value=0, step=BLUR_STEP_SIZE, title="Blur", width=100)
    multi_select = MultiSelect(title=output_filename, value=[default_term], options=options, width=100, height=300)
    checkbox = CheckboxGroup(labels=["TFIDF"], active=[])
    radio_button_group = RadioGroup(labels=["Viridis256", "Greys256"], active=0)
    button = Button(label="Metadata",button_type="default", width=100)
    stats = Button(label="Statistics",button_type="default", width=100)
    if class_id:
        label = "Show "+str(class_id)
        checkbox_class = CheckboxGroup(labels=["Show "+str(class_id)], active=[])

    # Javacode
    code_callback_slider1 = return_JS_code('slider1')
    code_callback_slider2 = return_JS_code('slider2')
    code_callback_ms = return_JS_code('multi_select')
    code_callback_checkbox = return_JS_code('checkbox')
    code_callback_rbg = return_JS_code('rbg')
    code_callback_button = return_JS_code('button')
    code_callback_stats = return_JS_code('stats')
    if class_id:
        code_callback_class = return_JS_code('class')

    # Callbacks
    callback_slider1 = CustomJS(args={'source': source, 'mapper': mapper, 'slider2': slider2, 'checkbox': checkbox}, code=code_callback_slider1)
    callback_slider2 = CustomJS(args={'source': source, 'mapper': mapper, 'slider1': slider1, 'checkbox': checkbox}, code=code_callback_slider2)
    callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source, 'slider1': slider1, 'slider2': slider2, 'checkbox': checkbox, 'p': p, 'mapper': mapper}, code=code_callback_ms)
    callback_checkbox = CustomJS(args={'source': source, 'slider2': slider2, 'mapper': mapper, 'slider1': slider1}, code=code_callback_checkbox)
    callback_radio_button_group = CustomJS(args={'p': p, 'multi_select': multi_select, 'mapper': mapper, 'term_to_class': term_to_class, 'Viridis256': Viridis256, 'Greys256': Greys256}, code=code_callback_rbg)
    callback_button = CustomJS(args={'term_to_metadata': term_to_metadata, 'multi_select': multi_select},code=code_callback_button)
    callback_stats = CustomJS(args={'term_to_stats': term_to_stats, 'multi_select': multi_select},code=code_callback_stats)
    if class_id:
        callback_ms = CustomJS(args={'source': source, 'term_to_source': term_to_source, 'slider1': slider1, 'slider2': slider2,
        'checkbox': checkbox, 'p': p, 'mapper': mapper, 'source_class': source_class, 'term_to_class': term_to_class, 'class_hex': class_hex,
        'checkbox_class': checkbox_class}, code=code_callback_ms)
        callback_class = CustomJS(args={'multi_select': multi_select, 'term_to_class': term_to_class, 'class_hex': class_hex}, code=code_callback_class)

    # On change
    slider1.js_on_change('value', callback_slider1)
    slider2.js_on_change('value', callback_slider2)
    multi_select.js_on_change("value", callback_ms)
    checkbox.js_on_change('active', callback_checkbox)
    radio_button_group.js_on_change('active', callback_radio_button_group)
    button.js_on_event(events.ButtonClick, callback_button)
    stats.js_on_event(events.ButtonClick, callback_stats)
    if class_id:
        checkbox_class.js_on_change('active', callback_class)

    # Layout
    if class_id:
        layout = row(multi_select, p, column(slider1, slider2, checkbox, checkbox_class, radio_button_group, button, stats))
    else:
        layout = row(multi_select, p, column(slider1, slider2, checkbox, radio_button_group, button, stats))
    show(layout)

def get_tables(files):
    '''
    This function recieves a list of files, imports these files and their corresponding metadata, and saves them in a dictionary.
    '''
    tables = dict()
    for file in files:
        table = import_table(file)
        term = file.split('/')[1].split('_')[0]
        metadata_file = 'metadata/'+str(term)+'.txt'
        metadata_lines = open(metadata_file, 'r')
        metadata = metadata_lines.readlines()

        tables[term] = {'table': table, 'metadata': metadata}

    return tables

def get_files(folder):
    '''
    This function recieves the input folder and returns a list of file paths in this folder.
    '''
    files = ['%s/%s' % (folder, file) for file in listdir(folder) if '.pkl' in file]
    return files

def parser():
    parser = argparse.ArgumentParser(description='This script makes a table of the query IDs, their names and their properties')
    parser.add_argument('-i', required=True, metavar='input_folder', dest='input_folder', help='[i] to select input folder')
    parser.add_argument('-o', required=True, metavar='output_filename', dest='output_filename', help='[o] to name .html output filename (a \'plots\' folder is required!)')
    parser.add_argument('-xmin', required=False, metavar='xmin', dest='xmin', help='[xmin] to select x axis minimum (logP), default is -5')
    parser.add_argument('-xmax', required=False, metavar='xmax', dest='xmax', help='[xmax] to select x axis maximum (logP), default is 10')
    parser.add_argument('-ymax', required=False, metavar='ymax', dest='ymax', help='[ymax] to select y axix maximum (mass in Da), default is 1600')
    parser.add_argument('-class', required=False, metavar='class_id', dest='class_id', help='[c] to select a class to be shown in the plot with a click on the class button')
    arguments = parser.parse_args()
    return arguments

def main():
    # startTime = datetime.now()

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

    plot(tables, output_filename, xmin, xmax, ymin, ymax, args.class_id)

    # print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
