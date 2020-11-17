# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 10:43:07 2020

@author: wteasley

The purpose of this script is to combine and process dataframes sourced from
excel spreadsheets of ICP-MS data that have been exported from Agilent's 
Masshunter software.

The first several steps of processing ICP-MS or other elemental data is performed
and a series of excel spreadsheets demonstrating the data transformations as well 
as a reporting/quantification limit chart and a separate spreadsheet of method 
blanks are returned.

Future versions will further automate QC validation and data processing.

"""

import os
import re
from datetime import date
from itertools import combinations
from itertools import chain
import numpy as np
import pandas as pd


def concat_data(file_names):
    """
    Sources dataframes from excel spreadsheets of ICPMS data and concatenates
    them.  
    
    If column labels are not exactly the same between dataframes, pandas
    automatically sorts the columns alphabetically.  Data sets will generally
    vary with respect to which analytes were acquired.
    
    Columns are sorted as strings ascending alphabetically.  So don't be confused 
    that "103 Rh" is sorted before "23 Na".  This is addressed later with the
    function "reorder_df".
    
    The algorithm below starting around line 50 (as of 11/10/20) is designed to 
    identify groups of dataframes that should be combined due to the 
    "calibration trick" used for acquisitions/elements with separate groups of 
    calibration standards.  This prevents duplicate data in the final dataframe.
    
    """
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
    frames = [pd.read_excel(workbook, header=[0,1]) for workbook in file_names]
    
    # # For debugging.
    # for i, frame in enumerate(frames):
    #     frame.to_excel(str(i) + ".xlsx")
    
    ###
    
    # The "calibration trick" (see above) needed for certain analytes makes this
    # algorithm necessary.

    # Create an iterator of dataframe combinations.
    # Note: iterator destroyed once iterated through or used with function.
    comb_df = combinations(frames, 2)
    # Create an iterator of index position combinations.
    idx_positions = []
    for i, v in enumerate(frames):
        num = int(i)
        idx_positions.append(num)
    comb_idx = combinations(idx_positions, 2)
    
    # Create an dict to group together dataframes to be combined, using the
    # first date-time value as the key and the index positions of the dataframes
    # within the list 'frames' as values.  Checks the date-time column of each
    # dataframe and groups the dataframes together if the check column is identical.
    dfs = {}
    for i, v in zip(list(comb_idx), list(comb_df)):
        # For debugging:
        # print(i)
        # print(f"{v[0].shape}, {v[1].shape}")
        chk_col1 = v[0][('Sample', 'Acq. Date-Time')]
        chk_col2 = v[1][('Sample', 'Acq. Date-Time')] 
        # Make a new dict entry for each df group to be combined.
        # Use pd.equals() to avoid ValueError using "==".
        if chk_col1.equals(chk_col2):
            if chk_col1[0] not in dfs.keys():
                dfs[chk_col1[0]] = set()
                dfs[chk_col1[0]].add(int(i[0]))
                dfs[chk_col1[0]].add(int(i[1]))
            else:
                dfs[chk_col1[0]].add(int(i[0]))
                dfs[chk_col1[0]].add(int(i[1]))
    
    # Combine each group of dataframes into a single dataframe.
    result = []            
    for k in dfs.keys():
        count = 0
        max_count = len(dfs[k]) - 1
        for i in dfs[k]:
            #print(f"Count: {count}")
            #print(f"Max Count: {max_count}")
            if count == 0:
                #print("First combine.")
                new = frames[i].combine_first(frames[i+1])
            elif i < max_count:
                #print("Next combine.")
                new = new.combine_first(frames[i+1])
            count += 1
        result.append(new)
    
    # Create final dataframe list for concatenating.
    # Can't use .remove() on dataframe objects in list.
    final = [df for i, df in enumerate(frames) \
             if i not in flatten_dict_values(dfs)]
    final = final + result
    
    ###
    
    return pd.concat(final, ignore_index=True)


def flatten_dict_values(dictionary):
    '''
    Subfunction.
    Used in function 'concat_data' to flatten dict values / list of sets.
    https://lothiraldan.github.io/2011-07-07-flatten-dictionaries-values/
    Doesn't work unless function is returned as list. 
    '''
    return list(chain(*dictionary.values()))


def reorder_df(df):
    """
    Reorder the columns after concatenation, which may cause the columns to be 
    undesirably reordered alphabetically.
    
    Uses top-level column label "Sample", all columns of which contain
    sample information, as a key for the placement of the first columns, followed
    by the elements in ascending mass order.
    """
    
    # Sort by Acq. Date-Time in case it got jumbled in function 'concat_data'.
    df = df.sort_values(('Sample', 'Acq. Date-Time'))
    
    # Put columns with "Sample" (sample ID information) header first.
    lst1 = [i for i in df.columns if i[0] == 'Sample']
    lst2 = [i for i in df.columns if i[0] != 'Sample']
    # Extract mass number from strings to put masses in proper order.
    lst3 = [int(i[0][:3]) for i in lst2]
    paired_lst = list(zip(lst3, lst2))
    paired_lst.sort()
    lst4 = [i[1] for i in paired_lst]
    order = lst1 + lst4
    
    return df[order]
    
    
def process_df(df):
    """
    The majority of data munging is carried out in this function.

    Users will likely need to tailor this section to thier own data, but the 
    general process should be relevant.

    There is an assumption that the analyst is verifying certain QC measures 
    during the run, so they will be removed at this point (e.g. check standards).

    """
    
    # Replace problematic values.
    # Use NaNs as opposed to any sort of string placeholder value.
    pattern = r"<0.0+"
    df.replace(to_replace=pattern, value=0, inplace=True, regex=True)
    df.replace(to_replace="OR", value=np.NaN, inplace=True)

    # Remove columns of internal standard data.        
    lst = [(first, second) for first, second in df.columns if "ISTD" in first]
    df.drop(labels=lst, axis=1, inplace=True)
    
    # Remove other unwanted columns.
    unwanted = ['Data File', 'Level', 'Rjct', 
                'Unnamed: 0_level_1', 'Vial Number']
    lst = [('Sample', name) for name in unwanted]
    df.drop(labels=lst, axis=1, inplace=True)
    
    # Remove unwanted rows.
    # Decided, for now, that I don't want ICV, CCVs, Check blanks, and cal stds.
    wanted = ['Sample', 'Spike Ref', 'DUP', 'Spike QC']
    df = df[df[('Sample', 'Type')].isin(wanted)]

    # Remove rows without any data.
    col_lst = [col for col in df.columns if col[0] != 'Sample']
    df = df.dropna(axis=0, how='all', subset=col_lst)
    
    # # For debugging:
    # # Ensure all data is numeric (aside from sample information).
    # col_lst = [i for i in df.columns if i[0] != 'Sample']
    # # Gets SettingWithCopyWarning
    # for col in col_lst:
    #     for i in df.index:
    #         df[col][i] = pd.to_numeric(df[col][i])
    # # Doesn't get warning.
    # for col in col_lst:
    #     for i in df.index:
    #         pd.to_numeric(df[col][i])
    
    df = df.reset_index()

    return df


def widdle_to_raw_conc(df):
    
    # Select only raw ICPMS concentrations.
    lst1 = [col for col in df.columns if col[0] == 'Sample']
    lst2 = [col for col in df.columns if col[1] == 'Meas. Conc. [ ppb ]']
    
    return df[lst1 + lst2]


def method_blank_table(df):
    """
    Pulls rows from data table with "Method Blank" in either "Sample Name" or 
    "Comment" columns and creates a new dataframe.
    
    Method and reagent blanks should be removed prior to filtering values per 
    calibration range limits.

    (11/16/20) This will be eventually updated to recognize additional terms 
    such as "reagent blank" and variable casing by means of regular expressions.

    """
    mb = df[df[('Sample', 'Sample Name')].isin(['Method Blank'])]
    mb.reset_index(drop=True, inplace=True)
    
    # Remove method blank rows.
    # DOESN'T WORK - can't use .str on dataframe object.
    # df['Sample'].str.contains('Method Blank')
    # Doesnt work - for same reason above:
    # df[~df.loc[:,[('Sample', 'Sample Name'), 
    #               ('Sample', 'Comment')]].str.contains('Method Blank')]
    # WORKS - Generate boolean mask to take a slice from.
    # Need parameter "na=False" if NaNs
    df = df[~df.loc[:,('Sample', 'Sample Name')].str.contains(
        'Method Blank', na=False)]
    df = df[~df.loc[:,('Sample', 'Comment')].str.contains(
        'Method Blank', na=False)]
    
    df.reset_index(drop=True, inplace=True)
    
    return (mb, df)

def quant_limit_table(df):
    '''
    Creates a quantification limit table based on the elements acquired in the 
    dataset and a csv file "quant_limit_table_ICPMS.csv" with commonly used 
    values.  Note that these are not true analytical limits of quantification, 
    but rather limits the analyst can use to filter and flag data based off of
    the calibration range limits and knowledge of the instruments linear 
    dynamic range.
    
    The csv file needs to be in same directory as the data sets.

    The function uses partial string matching so that a column label like 
    "78 Se [He gas]" is recognized as selenium and a corresponding column is 
    added to the quant limit table.  Also, if two isotopes of selenium are 
    acquired, the resulting quant limit chart will have two columns for selenium.
    
    Note that in the csv file, column titles have a space character on either 
    side to prevent undesired partial string matching (e.g. " Se ").  This may
    need to be altered per user data.

    '''    
    # Need to specify index column (string values).
    # For some reason, read_csv reads in blank columns.  Drop them.
    quant_limit = pd.read_csv("quant_limit_chart_ICPMS.csv", index_col=0)
    quant_limit.dropna(axis='columns', how='all', inplace=True)
    # Rename for convenience.
    ql = quant_limit
    
    ###
    # Select elements from quant limit chart based on elements in main df.
    # Some attempts to find partial matches of column labels:
    #df.columns.str.contains('Na')
    # AttributeError: Can only use .str accessor with Index, not MultiIndex
    #df.columns[5][0].str.contains('Na')
    # AttributeError: 'str' object has no attribute 'str'
    # Logic that can be used:
    #ql.columns[1] in df.columns[5][0]
    # True
    ###
    
    ###
    # # Works with first return "return ql[lst]"
    # # Create a simple list to avoid working with multi-index.
    # df_col = [col[0] for col in df.columns]
    # # Create a list of elements contained the dataset.
    # lst = [col2 for col1 in df_col for col2 in ql.columns if col2 in col1]
    
    # 2nd attempt: 
    # Giving same column label to matching column in ql table.
    # This is so dynamic matching can be used to work between tables later.
    # Create a list to avoid working with multi-index.
    df_col = [col[0] for col in df.columns]
    lst1 = [col2 for col1 in df_col for col2 in ql.columns if col2 in col1]
    # Create a list of elements contained the dataset.
    lst2 = [col1 for col1 in df.columns for col2 in ql.columns \
            if col2 in col1[0]]
    ###
    
    # Modify quant limit table so that columns (elements) match.
    ql1 = ql[lst1]
    
    # Rename columns in quant limit table to match target dataframe.
    # Can't use dict because of requirement to have unique keys (there can't -
    # be two instances of the same element, like if you acquired two isotopes).
    # Can't use .rename() because if 2 instances of the same column name - 
    # they will be renamed simultaneously with duplicate results.
    ql1.columns = lst2
    
    #return ql[lst]
    return ql1

def method_blank_qc(df, mb, ql):
    '''
    Not functioning as of 11/16/20.

    Should only be one method blank per date. (see comment at line)
    
    '''
    
    col = ql.columns[0]
    # Set the value to something different so you can test it.
    ql.at['Low - Flag', col] = 1
    # Only replaces value if condition is true (works though).
    mb.loc[:,col].mask(mb.loc[:,col] > ql.at['Low - Flag', col], True)
    
            
    ###
    
    # For testing.
    df = df_no_blank
    mb_test_values = [1, 0.1, 0.1, 1, 1, 0.1, 0.1, 1, 0.1, 0.1, 1, 0.1, 0.1, 
                      0.1, 0.1, 0.1, 0.1]
    ql.loc['Low - Flag'] = mb_test_values
    
    # Test if any method blank results are over corresponding low cal. limit.
    col_lst1 = [col for col in mb.columns if col[0] == 'Sample']
    col_lst2 = [col for col in mb.columns if col[0] != 'Sample']
    # Compare mb to ql column-wise; generate boolean mask.
    # Note that NaNs register as false.
    mb_bool = mb[col_lst2] > ql.loc['Low - Flag']
    
    ###
    # 1st approach:
    df[col_lst2].where(mb_bool, 'FLAG')
    
    # 2nd approach:
    mb_bool_whole = pd.concat([mb[col_lst1], mb_bool], axis=1)
    # ERROR "ValueError: Boolean array expected for the condition, not object"
    # Can't use a partial boolean array (sample information included).
    df.where(mb_bool_whole, 'FLAG')
    
    # 3rd approach:
    # Convert date time to date for index-matching and drop unwanted
    # sample information columns (that won't become multi-index).
    # BETTER WAY TO PARSE DATE AND TIME INTO SEPARATE COLUMNS?
    df = df_no_blank # FOR TESTING
    s = df.loc[:,('Sample','Acq. Date-Time')]
    s.apply(lambda s: s.date())
    df[('Sample', 'Date')] = s.apply(lambda s: s.date())
    df = df[col_lst1 + [('Sample', 'Date')] + col_lst2]
    drop = [('Sample', 'Acq. Date-Time'), ('Sample', 'Total Dil.'),
            ('Sample', 'Type'), ('Sample', 'Comment'), ('Sample', 'Sample Name')]
    df.drop(columns=drop, inplace=True)
    # Add date column to mb
    s = mb.loc[:,('Sample','Acq. Date-Time')]
    s.apply(lambda s: s.date())
    mb[('Sample', 'Date')] = s.apply(lambda s: s.date())
    mb = mb[col_lst1 + [('Sample', 'Date')] + col_lst2]
    drop = [('Sample', 'Acq. Date-Time'), ('Sample', 'Total Dil.'),
            ('Sample', 'Type'), ('Sample', 'Comment'), ('Sample', 'Sample Name')]
    mb.drop(columns=drop, inplace=True)
    # Create multi-index for columns.
    df = df.set_index([('Sample', 'Date')])
    mb = mb.set_index([('Sample', 'Date')])
    # Verify index values are the same.
    mb.index.unique() == df.index.unique()
    # Create boolean mask.
    mb_bool = mb[col_lst2] < ql.loc['Low - Flag']
    df.where(mb_bool, 'FLAG')
    # ValueError: operands could not be broadcast together with shapes (85,17) (49,17) () 
    # Need to use .mask() because you only want to replace values that -
    # correspond to "True" (and not False as .where() does)?
    
    # 4th Approach:
    # Need to make 'date' columns first.
    df_test = df
    mb_bool_whole = pd.concat([mb[col_lst1], mb_bool], axis=1)
    # Returns a datetime timestamp
    df_test[('Sample', 'Date')].iloc[0]
    # Returns a datetime date.
    df_test[('Sample', 'Acq. Date-Time')].iloc[0].date()
    # Returns only element columns from mb_bool table for 7/13/20 run.
    mb_bool_whole[col_lst2].loc[mb_bool_whole[
        ('Sample', 'Date')] == date(2020,7,13)]
    
    s = mb_bool_whole.loc[:,('Sample','Acq. Date-Time')]
    s.apply(lambda s: s.date())
    mb_bool_whole[('Sample', 'Date')] = s.apply(lambda s: s.date())
    mb_bool_whole = mb_bool_whole[col_lst1 + [('Sample', 'Date')] + col_lst2]
    df.columns == mb_bool_whole.columns
    # ERROR: IndexError: only integers, slices (`:`), ellipsis (`...`), 
    # numpy.newaxis (`None`) and integer or boolean arrays are valid indices.
    # For some reason, iloc returns a series while loc returns a dataframe.
    # Shape of loc return is (1, 17) - why would it return dataframe?
    for i in df_test.index:
        df_test[col_lst2].iloc[i] = mb_bool_whole[col_lst2].loc[
            mb_bool_whole[('Sample', 'Date')] == df_test[
                ('Sample', 'Date')].iloc[i]]
    # 5th approach
    df = df_no_blank # FOR TESTING
    # Create an empty dataframe instead and append to it.
    mb_bool_whole = pd.concat([mb[col_lst1], mb_bool], axis=1)
    # Create and add date column to df.
    df = add_date_column(df)
    # Create and add date column to mb.
    mb_bool_whole = add_date_column(mb_bool_whole)
    # Create empty dataframe and append rows of boolean data to mirror df.                              
    test = pd.DataFrame()    
    for i in df_test.index:
        date = df_test[('Sample', 'Date')].iloc[i]
        col = mb_bool_whole[col_lst2].loc[
                mb_bool_whole[('Sample', 'Date')] == date]
        test = test.append(col)
    test.reset_index(drop=True, inplace=True)
    df_mb_bool = pd.concat([df_test[col_lst1], test], axis=1)
    
    df_mb_qc = df[col_lst2].mask(df_mb_bool[col_lst2], 'FLAG')
    
    return pd.concat([df[col_lst1], df_mb_qc], axis=1)


def add_date_column(df):
    #nonlocal col_lst1, col_lst2 # Do you need this?
    s = df.loc[:,('Sample','Acq. Date-Time')]
    s.apply(lambda s: s.date())
    df[('Sample', 'Date')] = s.apply(lambda s: s.date())
    return df[col_lst1 + [('Sample', 'Date')] + col_lst2]


def filter_quant_limits(df, ql):
    """
    Method/reagent blanks should be separated before executing this function.
    
    CURRENTLY SET TO FILTER "FLAG" VALUE.  Not flagging only filtering values
    for now for simplicity.
    
    Okay if double instance of same mass and column order doesn't matter -
    function in ql chart will make a chart with matching columns (elements).
    
    Works based on ql dataframe having same columns as df.
    
    """
    
    # Make a copy to prevent SettingWithCopy warning.
    df = df.copy()
    
    # Specify non-sample information columns to be altered.
    col_lst = [i for i in df.columns if i[0] != 'Sample']
    
    ###
    # 1st attempt:
    # # WARNING: SettingWithCopy
    # for col1 in col_lst:
    #     for col2 in ql.columns:
    #         for i in df.index:
    #             if df[col1][i] < ql[col2][1] or df[col1][i] > ql[col2][2]:
    #                 df[col1][i] = np.NaN
    
    # ql changed so that column labels match target df (no need to have - 
    # nested for loop).
    # CURRENTLY SET TO FILTER "FLAG" VALUE
    # DOES THIS WORK WITHOUT FOR LOOP???
    for col in col_lst:
        df.loc[:,col] = df.loc[:,col].where(
            df.loc[:,col] > ql.at['Low - Flag', col], np.NaN)
        
    for col in col_lst:
        df.loc[:,col] = df.loc[:,col].where(
            df.loc[:,col] < ql.at['High - Flag', col], np.NaN)
        
    # Drop any rows with no data.
    col_lst = [col for col in df.columns if col[0] != 'Sample']
    df = df.dropna(axis=0, how='all', subset=col_lst)
    
    # # For debugging (checking datatype).  9/7/2020.
    # # Produces array of data types in df, showing that all values are float.
    # # For some reason, array is transposed (df.shape is 70, 17, turns to 17,70)
    # # Can convert np.array to pd.DataFrame and export to csv for easier inspection.
    # lst = [[type(val) for val in df[
    #     col_lst].iloc[:,i]] for i in range(len(col_lst))]
    # arr = np.array([[type(val) for val in df[
    #     col_lst].iloc[:,i]] for i in range(len(col_lst))])
    
    return df


def main():

    # Date variable for generating files.
    d = date.today().strftime("%Y%m%d")
    
    file_names = ["20200713_ICPMS_FoliarSprays.xlsx", 
              "20200716_ICPMS_FoliarSprays.xlsx",
              "20200716_ICPMS_FoliarSprays_MgOnly.xlsx"]
    
    # Work
    os.chdir(r"F:\ICPMS-FoliarGardenSprays-PythonExercise")
    # Home
    #os.chdir(r"D:\ICPMS-FoliarGardenSprays-PythonExercise")
    
    df = concat_data(file_names)
    df.to_excel(d + "_df_s1_concat.xlsx")
    
    df_reordered = reorder_df(df)
    df_reordered.to_excel(d + "_df_s2_reorder.xlsx")
    
    df_processed = process_df(df_reordered)
    df_processed.to_excel(d + "_df_s3_processed.xlsx")
    
    df_raw = widdle_to_raw_conc(df_processed)
    df_raw.to_excel(d + "_df_s4_raw.xlsx")
     
    (mb, df_no_blank) = method_blank_table(df_raw)
    
    # Produces a matching table of quantification limits (high and low 
    # calibration limits) for each element and isotope in the main data set.
    ql = quant_limit_table(df_no_blank)
    
    df_final = filter_quant_limits(df_no_blank, ql)
    df_final.to_excel(d + "_df_s5_filter.xlsx")
    
    ql.to_excel(d + "_ql.xlsx")
    mb.to_excel(d + "_mb.xlsx")
    
if __name__ == '__main__':
    main()
    
