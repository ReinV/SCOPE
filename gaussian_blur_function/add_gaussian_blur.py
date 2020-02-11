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
    # columns for blur
    blur_columns_count = int(blur_max / step_size)

    # add blur columns to dataframe
    blur_list = []

    # prepare column names for the new dataframe
    column_names = []

    for column_name in df.columns:
        column_names.append(column_name)

    for i in range(0, blur_columns_count+1):
        sd_x = i*step_size
        if sd_x%1 == 0:
            sd_x = int(sd_x) # this to get the names of the column names equal to the bokeh slider values, e.g. slider returns 0 or 1 when python returns 0.0 or 1.0
        blur_list.append(sd_x)
        column_name = str(sd_x)
        column_names.append(column_name)

    # construct kernel with distances
    kernel_distances = {(-5,2):(7.5, sqrt(3)/2),(-5,3):(7.5, sqrt(3)/2),
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

    # make kernel for blur
    kernel_blur = {}
    for sd_x in blur_list:
        sd_y = sd_x / 2

        # calculate blur for all the keys in the kernel, replace distance with blur factor
        for key in kernel_distances:
            distance = kernel_distances [key]
            x = distance[0]
            y = distance[1]
            if sd_x == 0:
                blur = 0
            else:
                blur = get_blur(x, y, sd_x, sd_y)
            try:
                kernel_blur[key][sd_x] = blur
            except:
                kernel_blur[key] = {sd_x: blur}

    # Remake df in python dictionary
    coords_to_count = dict() # keeps all original coordinates and count values
    coords_to_blur = dict() # will have new coordinates and new count (blur) values

    # go through all coordinates in df
    for index, row in df.iterrows():
        q = row['q']
        r = row['r']
        counts = row['counts']
        tfidf = row['tfidf']
        coords_to_blur[(q,r)] = {}
        coords_to_count[(q,r)] = {'counts': counts, 'tfidf': tfidf}

        # add every column from the dataframe to the dictionary
        for i in range(len(df.columns)):
            column_name = df.columns[i]
            column_value = row[column_name]
            coords_to_blur[(q,r)][column_name] = column_value

    # put original count values in sd_x column
    for sd_x in blur_list:
        column_name = str(sd_x)
        # Add blur
        for coord in coords_to_count.keys():
            q = coord[0]
            r = coord[1]
            counts = coords_to_count[coord]['counts']
            tfidf = coords_to_count[coord]['tfidf']
            coords_to_blur[(q, r)][column_name] = [counts, tfidf]


    # do blur for every blur_column
    for sd_x in blur_list:
        column_name = str(sd_x)
        # Add blur
        for coord in coords_to_count.keys():
            q = coord[0]
            r = coord[1]
            counts = coords_to_count[coord]['counts']
            tfidf = coords_to_count[coord]['tfidf']

            # For every coordinate, apply kernel (to surrounding hexagons)
            for key in kernel_blur.keys():
                q_new = q + key[0]
                r_new = r + key[1]
                factor = kernel_blur[key][sd_x]
                counts_blur = counts*factor
                tfidf_blur = tfidf*factor

                # Check if surrounding hexagon exists and apply blur
                try:
                    coords_to_blur[(q_new, r_new)][column_name][0] += counts_blur
                    coords_to_blur[(q_new, r_new)][column_name][1] += tfidf_blur
                except:
                    # If hexagon does not exist, make a new hexagon
                    # All columns (see keys) need to be added here or else Bokeh will crash, but most values will be left empty
                    if sd_x == 0:
                        # make the row
                        coords_to_blur[(q_new, r_new)] = {'q': q_new, 'r': r_new, column_name: [counts_blur, tfidf_blur], 'counts': 0, 'tfidf': 0, 'scaling': 0, 'ids': [],
                        'id1': '', 'count1': '', 'name1': '',
                        'id2': '', 'count2': '', 'name2': '',
                        'id3': '', 'count3': '', 'name3': '',
                        'id1_tfidf': '', 'count1_tfidf': '', 'name1_tfidf': '',
                        'id2_tfidf': '', 'count2_tfidf': '', 'name2_tfidf': '',
                        'id3_tfidf': '', 'count3_tfidf': '', 'name3_tfidf': ''}
                    else:
                        # row should already be there, only need to add column value
                        coords_to_blur[(q_new, r_new)][column_name] = [counts_blur, tfidf_blur]
                        # coords_to_blur[(q_new, r_new)][column_name] = [counts_blur, tfidf_blur]


    df_gaussian = pd.DataFrame()

    # return dictionary to pandas dataframe
    for column_name in column_names:
        # add all dictionary keys as columns (lists) to the dataframe
        column_list = []
        for coords in coords_to_blur.keys():
            column_value = coords_to_blur[coords][column_name]
            column_list.append(column_value)
        # add column to the df
        df_gaussian[column_name] = column_list

    df_gaussian = df_gaussian.drop("ids", axis=1)
    return df_gaussian
