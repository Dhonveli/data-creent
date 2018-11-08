#!/usr/bin/env python

import os
import sys
import argparse
import shutil
import csv
import pandas as pd
import gzip

# Parses the command line arguments of this program


def parse_args():

    # Creates an object to parse the command line parameters
    parser = argparse.ArgumentParser(description="This program maps genes name to nessra code",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Positional parameters

    # Optional parameters
    parser.add_argument("-i", "--input", nargs=1,
                        help="Folder containing subfolder with expression profiles downloaded from TCGA", type=str)
    parser.add_argument("-a", "--annotation", nargs = 1, 
                        help="Annotation file to filter genes",type=str,default="hgnc_cc_zero_filtered_anno.csv")
    # Parses and returns the object containing the params
    return parser.parse_args()

def load_files(args):

    global array_matrix

    n = -1
    mulfile = []
    anno = []
    for root,dirs,files in os.walk(args.input[0]):
        n += 1
        if n == 0:
            dir_list = dirs
        if n < 3:
            for el in files:
                if el.endswith(".gz"):
                    sinfile = pd.read_csv(args.input[0]+"/"+dir_list[n-1]+"/"+el,
                    compression="gzip",sep="\t",index_col=0,header=None,names=[el.split(".")[0]])
                    mulfile.append(sinfile)
    array_matrix = pd.concat(mulfile,axis=1)
    array_matrix.to_csv("TCGA-PRAD", sep=',',header=1)
                    

def filter_matrix()


if __name__ == "__main__":
    
    # Parser
    args = parse_args()

    load_files(args)