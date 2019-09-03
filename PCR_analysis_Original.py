#!/usr/bin/env python

# Author: Eric Helfrich
# Based off of the original code of Abrar Al-Shaer
# This program filters unwanted PCR sample results and averages technical replicates within the file from
# the mean and standard deviation (blank cells were removed in excel).
# After filteration the program calculates the Ct(also known as Cq) value for each sample normalized to the reference,
# fold change, SEMs, and averaged standard deviations.

import pandas as pd
from pandas import DataFrame
import numpy as np
from scipy import stats
import os
import sys
import time

from PCR_functions import *


def main():
    reference = input("What is your reference?: ")
    print("Your Reference: ", reference) #  GAPDH
    one_target = input("What is your target?: ")
    print("Your target: ", one_target) #  AvUCP

    # Create list of data files and check if correct
    data_files = []
    directory = "input"
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".csv"):
                data_files.append(file)

    print(data_files)
    if yes_no("Are the data files listed above correct? [y/n]: "):
        print("Using printed files")
    else:
        print("Please correct files... \n.\n.\n...Exiting Program")

    # Check to see if merge process needs to be run
    merge = False
    if yes_no("Would you like to calculate ddCt for these datasets? [y/n]: "):
        print("ddCt calculations will be run")
        merge = True
        time.sleep(3)
    else:
        print("ddCt calculations will NOT be run")
        time.sleep(3)

    # Get the Cq Calculations csv file name
#    file_header = input("What do you want the Cq_calculations output files to start with?\n"
#                        "Format: [fileName]_Cq_calculations_[fileNum].csv")

    # Call functions for all PCR files
    sorted_dfs = []  # declare empty list for the results of the Ct_calculations_print function
    i = 0  # index for number of file
    for data in data_files:
        print(data)
        i = i + 1
        csv = csv_init(os.path.join(directory, data))
        target, content, CqAvg, CqDev, sample, set_name = rows_init_store(csv)
        sampleAll, set_nameAll, CqAvgAll, CqDevAll = Ct_calculations(target, content, CqAvg, CqDev,
                                                                     sample, set_name, reference, one_target)
        # Prints the results to screen and file, also creates a list of DFs for the merge process
        calcs = Ct_calculations_print(sampleAll, set_nameAll,
                              CqAvgAll, CqDevAll, str(i), data[:-4])
        sorted_dfs.append(calcs)

    #  ddCt Calculations process

    if merge:  # Only runs if merge is True

        input_cal = input("What is your calibrator?: ")

        # Calling function for merger of PCR Files and mean/SEM calculation for each set
        print("\n********PCR Calculations Merge*************\n")
        df = means_sem_calculation(sorted_dfs)
        print("\n********PCR Averages & SEM*************\n")
        print(df)

        # Calling function for calculations delta delta Ct
        print("\n********PCR Delta Delta Ct*************\n")
        df = delta_delta_ct(df, calibrator=input_cal) #FpreA
        print(df)

        # Calling function for calculating fold change
        print("\n********PCR Fold Change*************\n")
        df = fold_change(df)

        # Final DataFrame
        print("\n********Fold Change File*************\n")
        print(df)  #print first 5 rows of dataframe
        df.to_csv(os.path.join("output", "fold_changes_output.csv"))

if __name__ == "__main__":
    main()
