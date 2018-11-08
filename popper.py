#!/usr/bin/env python

import os
import sys
import argparse
import shutil
import csv
import pandas as pd
import gzip
import mygene


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
        for el in files:
            if el.endswith(".gz"):
                sinfile = pd.read_csv(args.input[0]+"/"+dir_list[n-1]+"/"+el,
                compression="gzip",sep="\t",index_col=0,header=None,names=[el.split(".")[0]])
                mulfile.append(sinfile)
    array_matrix = pd.concat(mulfile,axis=1)
    array_matrix.index = [v.split('.')[0] for v in array_matrix.index]
    array_matrix.to_csv("TCGA-PRAD", sep=',',header=1)
                            

def filter_matrix():

    map_gene_to_anno = {}

    mg = mygene.MyGeneInfo()
    ens_anno = [id.split(".")[0] for id in array_matrix.index.values]
    symb_anno = mg.querymany(ens_anno , scopes='ensembl.gene', species='human')    
    symb_anno = pd.DataFrame(symb_anno)
    symb_anno = symb_anno[symb_anno['notfound'] != True]
    with open(args.annotation) as file:

        # Discard the first line
        file.readline()

        # Legge il csv (il modulo csv è necessario perché i dati
        # sono incapsulati nelle double quotes e sarebbe difficile leggerli altrimenti)
        for line in csv.reader(file, quotechar='"', delimiter=',',
                               quoting=csv.QUOTE_ALL, skipinitialspace=True):

            # Estrae il nome del gene e il codice dell'isoforma
            name = line[10]
            code = line[0]

            # Controlla dentro la mappa con il gene non ci sia ancora,
            # in tal caso genera una lista associata
            if name not in map_gene_to_anno:
                map_gene_to_anno[name] = []

            # Aggiunge questa riga nella mappa
            map_gene_to_anno[name].append(line)
    symb_anno_filt = symb_anno[symb_anno['symbol'].isin(map_gene_to_anno)]
    array_matrix_filt = pd.merge(array_matrix,pd.DataFrame(symb_anno_filt[["query","symbol"]]),left_index=True,right_on="query",how="inner")
    array_matrix_filt = array_matrix_filt.set_index('symbol')
    array_matrix_filt = array_matrix_filt.drop('query',1)
    array_matrix_filt.to_csv("TCGA-PRAD-filt", sep=',',header=1)

if __name__ == "__main__":
    
    # Parser
    args = parse_args()

    load_files(args)

    filter_matrix()