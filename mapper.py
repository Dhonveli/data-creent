#!/usr/bin/env python

import os
import sys
import argparse
import shutil
import csv

# Parses the command line arguments of this program
def parse_args ():

    # Creates an object to parse the command line parameters
    parser = argparse.ArgumentParser(description="This program maps genes name to nessra code",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Positional parameters
    

    # Optional parameters
    parser.add_argument("-i", "--input", nargs=1, help="input folder", type=str, default="networks")
    parser.add_argument("-a", "--annotation", help="annotation file from LBM1819", type=str, default="hgnc_cc_zero_filtered_anno.csv")
    parser.add_argument("-r", "--csvboinc", help="csv file from https://gene.disi.unitn.it/test/gene_h.php", type=str, default="genehome PC-IM history.csv")

    # Parses and returns the object containing the params
    return parser.parse_args()

# Mappa geni to code
map_gene_to_anno = {}
map_code_to_anno = {}
map_code_to_args = {}
header_runned = []

def load_files (args):

    global map_gene_to_anno
    global map_code_to_anno
    global map_code_to_args
    global header_runned
    global list_gene

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

            # Controlla dentro la mappa con il gene non c'è ancora,
            # in tal caso genera una lista associata
            if name not in map_gene_to_anno:
                map_gene_to_anno[name] = []

            # Aggiunge questa riga nella mappa
            map_gene_to_anno[name].append (line)

            # Controlla se non c'è il codice trovato
            # in tal caso genera una lista associata
            if code not in map_code_to_anno:
                map_code_to_anno[code] = []

            # Aggiunge questa riga nella mappa
            map_code_to_anno[code].append (line)

    # Carica il file dei geni già runnati su boinc
    with open(args.csvboinc) as file:

        # Parsa come csv
        csvfile = csv.reader(file, quotechar='"', delimiter=',',
                    quoting=csv.QUOTE_ALL, skipinitialspace=True)

        # Preleva l'header che lo usera nel file already-runned.csv
        header_runned = next(csvfile)

        # Scorre linea per linea
        for line in csvfile:
            
            # Prende solo i geni i file che iniziano con 'T'
            # perché significa che sono singol gene expansion
            if line[2].startswith("T"):

                # Prende il codice dell'isoforma
                code = line[2].split('-')[0]

                # Inserisce nella mappa

                if code not in map_code_to_args:
                    map_code_to_args[code] = []
                
                map_code_to_args[code].append(line)

    

    
def process_input (gene_list, output_folder):

    global map_gene_to_anno
    global map_code_to_anno
    global map_code_to_args
    global header_runned
    global list_gene

    # Normalizza la path 
    output_folder = os.path.normpath(output_folder)

    try:
        # Rimuove la cartella di output
        shutil.rmtree (output_folder)

    except FileNotFoundError:
        pass

    # Carica la lista dei geni
    list_gene = [line.strip() for line in open(gene_list) if line.strip()]

    # Crea una nuova cartella
    os.mkdir (output_folder)

    # Nome della cartella onegene
    onegene = "onegenes"

    # Crea la cartella delle espansioni singoli genes
    os.mkdir (output_folder + '/' + onegene)

    # Scrive la local gene network
    with open(output_folder + '/hs-' + os.path.basename(gene_list) + '.txt', 'w') as filelgn:
        with open(output_folder + '/already-runned.csv', 'w') as filerunned:

            csvrunned = csv.writer(filerunned, quotechar='"', delimiter=',',
                     quoting=csv.QUOTE_ALL)

            csvrunned.writerow(header_runned)

            filelgn.write('from,to')

            for gene in list_gene:

                if gene not in map_gene_to_anno:
                    print ("Gene " + gene + " not found in the annotation file.")

                else:
                    for isoform in map_gene_to_anno[gene]:

                        code = isoform[0]

                        # scrive questo codice nella lgn
                        filelgn.write('\n'+code+','+code)

                        # Decide il nome dei file onegene
                        filename = code + '-' + gene + '.txt'

                        # controlla se è già stato runnato
                        if code in map_code_to_args:
                            for row in map_code_to_args[code]:
                                csvrunned.writerow (row)

                        # Altrimenti scrive il file da mandare a walter
                        else:
                            with open(output_folder + '/' + onegene + '/' + filename, 'w') as file:
                                file.write('from,to\n'+code+','+code)

            




if __name__ == "__main__":

    # Parses command line parameters
    args = parse_args()

    load_files(args)

    # Normalizza la path 
    args.input = os.path.normpath(args.input)

    for name in os.listdir(args.input):
        
        filename = args.input + '/' + name

        if os.path.isfile(filename):
            process_input(filename, filename + '_out')