#!/usr/bin/env python

import os
import sys
import argparse
import shutil
import csv

# Parses the command line arguments of this program


def parse_args():

    # Creates an object to parse the command line parameters
    parser = argparse.ArgumentParser(description="This program maps genes name to nessra code",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Positional parameters

    # Optional parameters
    parser.add_argument("-lgn", "--LocalGeneNetwork", nargs=1,
                        help="Folder containing LGN files", type=str)
    parser.add_argument("-sgn", "--SingleGeneExpansionList", nargs=1,
                        help="csv file containing gene list downloaded from OpenTarget", type=str)
    parser.add_argument("-a", "--annotation", help="annotation file from LBM1819",
                        type=str, default="hgnc_cc_zero_filtered_anno.csv")
    parser.add_argument("-r", "--csvboinc", help="csv file from https://gene.disi.unitn.it/test/gene_h.php",
                        type=str, default="genehome PC-IM history.csv")
    parser.add_argument(
        "-n", "--numbgenes", help="number of transcripts to run on NESSRA", type=int, default=50)

    # Parses and returns the object containing the params
    return parser.parse_args()

def load_files(args):

    global map_gene_to_anno
    global map_code_to_anno
    global map_code_to_args

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

            # Controlla se non c'è il codice trovato
            # in tal caso genera una lista associata
            if code not in map_code_to_anno:
                map_code_to_anno[code] = []

            # Aggiunge questa riga nella mappa
            map_code_to_anno[code].append(line)

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
            # perché significa che sono single gene expansion
            if line[2].startswith("T"):

                # Prende il codice dell'isoforma
                code = line[2].split('-')[0]

                # Inserisce nella mappa

                if code not in map_code_to_args:
                    map_code_to_args[code] = []

                map_code_to_args[code].append(line)

def process_lgn(gene_list):
    
    global map_code_to_anno
    global map_gene_to_anno
    global map_code_to_args
    global csvsummary
    global csvrunned

    list_gene = [line.strip() for line in open(gene_list) if line.strip()]

    with open(output_folder + '/hs-' + os.path.basename(gene_list) + '.txt', 'w') as filelgn:
        filelgn.write('from,to')
        # apre il file summary lgn
        # itera sui gene
        for gene in list_gene:
            # guarda se sono annotati
            if gene not in map_gene_to_anno:
                csvsummary.writerow(
                    [None, gene, "lgn", "not annotated",""])
            else:
                for isoform in map_gene_to_anno[gene]:

                    code = isoform[0]
                    # scrive questo codice nella lgn
                    filelgn.write('\n' + code + ',' + code)

                    # Decide il nome dei file onegene
                    filename = code + '-' + gene + '.txt'

                    # controlla se è già stato runnato
                    if code in map_code_to_args:
                        for row in map_code_to_args[code]:
                            csvrunned.writerow(row)
                            csvsummary.writerow(
                                [code, gene, "lgn", "already run",""])

                    # Altrimenti scrive il file da mandare a walter
                    else:
                        csvsummary.writerow(
                            [code, gene, "lgn", "to be run","HIGH"])
                        with open(output_folder + '/' + filename, 'w') as file:
                            file.write('from,to\n' +
                                        code + ',' + code)


def process_sgn(gene_list):
    global map_code_to_anno
    global map_gene_to_anno
    global map_code_to_args
    global csvsummary
    global csvrunned
    global numbgenes
    gene_list.readline()
    gene_rank = {}

    for line in gene_list:
        tmp = line.split(",")
        gene_rank[tmp[0]] = (float(tmp[3]) + float(tmp[4])) / 2
    sorted_gene = sorted(gene_rank, key=gene_rank.get, reverse=True)
    
    for gene in sorted_gene:
        if numbgenes <= 0:
            break
        else:
            if gene not in map_gene_to_anno:
                csvsummary.writerow(
                    [None, gene, str(gene_rank[gene]), "not annotated",""])
            else:
                for isoform in map_gene_to_anno[gene]:

                    code = isoform[0]

                    # Decide il nome dei file onegene
                    filename = code + '-' + gene + '.txt'

                    # controlla se è già stato runnato
                    if code in map_code_to_args:
                        for row in map_code_to_args[code]:
                            csvrunned.writerow(row)
                            csvsummary.writerow(
                                [code, gene, str(gene_rank[gene]), "already run",""])

                    # Altrimenti scrive il file da mandare a walter
                    else:
                        numbgenes -= 1
                        with open(output_folder + '/' + filename, 'w') as file:
                                file.write('from,to\n' +
                                            code + ',' + code)
                        if gene_rank[gene] > 0.6:
                            csvsummary.writerow(
                                [code, gene, str(gene_rank[gene]), "to be run","HIGH"])
                        else:
                            csvsummary.writerow(
                                [code, gene, str(gene_rank[gene]), "to be run","LOW"])
   

# Mappa geni to code
map_gene_to_anno = {}
map_code_to_anno = {}
map_code_to_args = {}

if __name__ == "__main__":

    # Parses command line parameters
    args = parse_args()

    load_files(args)

    numbgenes = args.numbgenes

    input = False
    if os.path.exists("output_folder"):
        shutil.rmtree("output_folder")

    # Normalizza la path
    os.mkdir("output_folder")

    # Normalizza la path
    output_folder = os.path.normpath(
        'output_folder')  # removes redundant separator
    
    # apre summary
    with open(output_folder + '/summary.csv', "w") as summary:
        # scrive header
        csvsummary = csv.writer(summary, quotechar='"', delimiter=',',
                                    quoting=csv.QUOTE_ALL)
        csvsummary.writerow(["ID", "GENE", "SCORE", "STATUS","PRIORITY"])

        # apre already_runned
        with open(output_folder + '/already-runned.csv', 'w') as filerunned:
            csvrunned = csv.writer(filerunned, quotechar='"', delimiter=',',
                                        quoting=csv.QUOTE_ALL)
            csvrunned.writerow(["ID","ORG","LGN","EXP","LastUpd","alpha","tsize","iter","nPC","nWU"])
            
            if args.LocalGeneNetwork is not None:
                input = True
                args.LocalGeneNetwork = os.path.normpath(args.LocalGeneNetwork[0])

                for name in os.listdir(args.LocalGeneNetwork):  # estrae i vari lgn

                    filename = args.LocalGeneNetwork + '/' + name

                    if os.path.isfile(filename):
                        process_lgn(filename)  # processo lgn

            if args.SingleGeneExpansionList is not None:
                input = True
                args.SingleGeneExpansionList = os.path.normpath(
                    args.SingleGeneExpansionList[0])

                with open(args.SingleGeneExpansionList, "r") as filesg:
                    process_sgn(filesg)
            else:
                if input == False:
                    raise Exception("No input was given!")

