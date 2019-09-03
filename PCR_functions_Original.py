# Author: Eric Helfrich
# Based off of the original code of Abrar Al-Shaer
# This file includes the the functions originally created by Abrar and updated to python 3
# and changed to meet updated requirements needed by the Smith Lab and UNC

# Imports
import pandas as pd
from scipy import stats
import os


def csv_init(file_PCR):
    print("---------Initiation of CSV file & Filteration---------")
    columns = []  # flag = columns
    from_csv = pd.read_csv(file_PCR, sep=None, header=0, engine='python')  # reading in the CSV file
    columns = from_csv[
        ['Target', 'Content', 'Sample', 'Biological Set Name', 'Cq', 'Cq Mean', 'Cq Std. Dev']]  # reads header columns
    # filteration step
    # blank Cq rows (or NTC) were removed in excel via F5 command
    # df = dataframe
    df = columns[columns.Cq != 0]  # removes zeroes if there is a zero in the Cq column
    df = df.drop_duplicates('Content', keep='first')  # removing duplicates based on content column value
    print(df)  # print first 5 lines of dataframe
    return df  # function returns dataframe as output


def rows_init_store(dataframe):
    print("---------Defining and Storing Rows From Dataframe---------")
    # Initializing blank lists
    Cq = []
    CqDev = []
    CqAvg = []
    content = []
    target = []
    sample = []
    set_name = []
    # loop through all the rows in the dataframe and add each elemenof each row to the proper column list
    for index, row in dataframe.iterrows():
        print(row[1], row[0], row[5], row[6])  # print content variable, target, mean, and standard deviation
        target.append(row[0])  # append = add to list
        content.append(row[1])
        Cq.append(row[4])
        CqAvg.append(row[5])
        CqDev.append(row[6])
        sample.append(row[2])
        set_name.append(row[3])
    return target, content, CqAvg, CqDev, sample, set_name  # function returns lists as output


def Ct_calculations(target, content, CqAvg, CqDev, sample, set_name, reference, one_target):
    print("--------------Ct Calculation Inputs---------------")
    CqAvgAll = []
    CqDevAll = []
    sampleAll = []
    set_nameAll = []
    reference_flag = 0
    target_flag = 0
    if reference in target[0]:
        reference_flag = 1
    if one_target in target[0]:
        target_flag = 1
    print("Reference FLAG:", reference_flag)
    print("Target FLAG:", target_flag)
    for i in range(0, len(content), 2):
        if reference_flag == 1:
            print(sample[i], target[i + 1], "-", target[i], "Means:", CqAvg[i + 1], "-", CqAvg[i],
                  "Standard Dev:", CqDev[i + 1], "&", CqDev[i])
            sampleAll.append(sample[i])
            set_nameAll.append(set_name[i])
            CqAvgAll.append(float(CqAvg[i + 1]) - float(CqAvg[i]))  # Cq calculation wither averages
            CqDevAll.append(
                (float(CqDev[i + 1]) ** 2.0 + float(CqDev[i]) ** 2.0) ** 0.5)  # Cq calculation with standard deviations
        if target_flag == 1:
            print(sample[i], target[i], "-", target[i + 1], "Means:", CqAvg[i], "-", CqAvg[i + 1],
                  "Standard Dev:", CqDev[i], "&", CqDev[i + 1])
            sampleAll.append(sample[i])
            set_nameAll.append(set_name[i])
            CqAvgAll.append(float(CqAvg[i]) - float(CqAvg[i + 1]))  # Cq calculation wither averages
            CqDevAll.append(
                (float(CqDev[i]) ** 2.0 + float(CqDev[i + 1]) ** 2.0) ** 0.5)  # Cq calculation with standard deviations
    return sampleAll, set_nameAll, CqAvgAll, CqDevAll


def Ct_calculations_print(sampleAll, set_nameAll, CqAvgAll, CqDevAll, fileNum, header):
    print("------------Ct Calculation Results-----------------")
    for i in range(0, len(CqAvgAll)):  # loop through and print all the contents of the lists
        print(sampleAll[i], set_nameAll[i], CqAvgAll[i], CqDevAll[i])
    # Ct Calculations DataFrame
    print("\n********Ct Calc File*************\n")
    df_Ct = pd.DataFrame({'Sample IDs': sampleAll, 'Biological Sets': set_nameAll, 'Cq Averages': CqAvgAll,
                          'Cq Standard Deviations': CqDevAll})  # creates a dataframe
    print(df_Ct)  # print all rows of dataframe
    print("\n********Ct Calc File SORTED*************\n")
    df_sort = df_Ct.sort_values('Biological Sets')
    print(df_sort)
    #  make a CSV file out of the DataFrame unique to the file number (fileNum)
    file_name = header + "_NUM_" + fileNum + "_Cq_calculations_" + '.csv'
    df_Ct.to_csv(os.path.join("output", file_name))  # print CSV contents to a file
    return df_sort


def means_sem_calculation(df_sorted):
    """
    Merge all the files and calculate the average and SEM for each biological set.
    :param df_sorted: List of DataFrame
    :return: DataFrame
    """
    merged_df = pd.DataFrame(columns=['Sample IDs', 'Biological Sets', 'Cq Averages', 'Cq Standard Deviations'])
    for df in df_sorted:
        merged_df = merged_df.append(df, ignore_index=True, sort=False)

    means = merged_df.groupby('Biological Sets')['Cq Averages'].mean()
    # Group by the values in Biological Sets and apply the sem function to the values in the group
    SEMs = merged_df.groupby('Biological Sets')['Cq Averages'].apply(lambda x: stats.sem(x.values))

    final_df = pd.DataFrame(dict(Mean = means, SEM = SEMs)).reset_index()

    return final_df


def delta_delta_ct(input_df, calibrator):
    """
    Calculate the delta delta ct for each biological set.
    :param input_df: DataFrame with the mean/SEM of each biological set
    :param calibrator: What variable that is the calibrator for the dataset.  Set by user at runtime
    :return: DataFrame
    """

    calibrator_mean = input_df['Mean'].loc[input_df['Biological Sets'] == calibrator].values[0]

    input_df['ddCt'] = input_df['Mean'].apply(lambda x: x - calibrator_mean)

    return input_df


def fold_change(input_df):

    input_df['FC'] = input_df['ddCt'].apply(lambda x: 2 ** (x * -1))
    input_df['FC Upper'] = input_df.apply(lambda x: 2 ** ((x['ddCt'] + x['SEM']) * -1), axis=1)
    input_df['FC lower'] = input_df.apply(lambda x: 2 ** ((x['ddCt'] - x['SEM']) * -1), axis=1)

    return input_df

def yes_no(question):
    response = input(question + "(y/n): ").lower().strip()
    print("")
    while not (response == "y" or response == "yes" or response == "n" or response == "no"):
        print("Input yes or no")
        response = input(question + "(y/n):").lower().strip()
        print("")
    if response[0:] == "y" or response[0:] == "yes":
        return True
    else:
        return False

