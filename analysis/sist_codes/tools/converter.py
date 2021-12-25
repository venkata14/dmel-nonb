#!/usr/bin/python
import argparse
import numpy as np
import pandas as pd
import glob
import os

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputdir", help="The input sequence dir for one temperature")
parser.add_argument("-o", "--outputfile", help="This is .csv outputfile for the merged files")
args = parser.parse_args()



def parse_sist(fn):
    """Parse raw SIST output; return an array of positions and probabilities"""
    data = []
    with open(fn,'r') as f:
        for line in f:
            line = line.strip()
            if 'Position' in line or 'WARNING' in line:
                continue
            line = line.split()
            # line = float(line)
            data.append(line)
    return np.array(data, dtype="float64")

def extractDF(window):
    """This takes in a WINDOW list (defined below) and outputs a dataframe of all the sist values in order. Remember that WINDOWS must be sorted first"""
    df = pd.DataFrame(columns=["Position", "P_melt", "P_Z", "P_cruciform"])
    for n in window:
        data = parse_sist(n)
        if data.shape[0] == 0:
            continue
        inner_df = pd.DataFrame(data, columns=["Position", "P_melt", "P_Z", "P_cruciform"])
        df = df.append(inner_df, ignore_index=True)
    return df




# Get all files from the inputdir
fns = glob.glob(
    os.path.join(args.inputdir, '*.sist')
)

# The files contain both the windows so have to seperate them
WINDOW1 = []
WINDOW2 = []

# Seperates the two windows into the two lists
for fn in fns:
    if "WINDOW1" in fn:
        WINDOW1.append(fn)
    elif "WINDOW2" in fn:
        WINDOW2.append(fn)


# Lazy Man's Bubble Sort to sort the files
for m in [WINDOW1, WINDOW2]:
    # ContIter = True
    L = len(m)
    for _ in range(L):
        for n in range(L):

            # Making sure there is no out of index error
            if n + 1 == L:
                continue

            # Checking if one file starts eariler or not and switching if they do
            first = int(m[n].split('-')[3])
            second = int(m[n+1].split('-')[3])
            if first > second:
                tmp = m[n]
                m[n] = m[n+1]
                m[n+1] = tmp
                # Method to turn off while loop
                # ContIter = True
            # else:
                # ContIter = False

# At this point both the WINDOW1 and WINDOW2 files are properly sorted and use extractDF

df1 = extractDF(WINDOW1)
df2 = extractDF(WINDOW2)

# Difference in length
dL = len(df1) - len(df2)
# Make df1 the bigger one so that I can add 0s to the smaller one (df2)
if dL < 0:
    tmp = df1
    df1 = df2
    df2 = tmp
    raise ValueError("df1 is supposed to be the bigger one. Did something go wrong?")

# Padding tho top of df2 so that it is the same length as df1
buffer = np.zeros((dL, 4))
bufferDF = pd.DataFrame(buffer, columns=["Position", "P_melt", "P_Z", "P_cruciform"])
df2  = bufferDF.append(df2, ignore_index=True)

# Next step is to create one dataframe with only the largest values from each index

# merging the two dataframes based on size
a = df1.values
b = df2.values
finalDF = pd.DataFrame(np.where(a > b, a, b), index=df1.index, columns=df1.columns)

# saving the dataframe to the output file name
print(args.outputfile)
finalDF.to_csv(args.outputfile, index=False)
