import os
import argparse


def parse_args():

    # Creates an object to parse the command line parameters
    parser = argparse.ArgumentParser(description="This program download resulting data from Gdrive and produce the txt files for the gene to be runned",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Positional parameters

    # Optional parameters
    parser.add_argument("-lgn", "--LocalGeneNetwork", nargs=1,
                        help="Folder containing LGN files",
                        type=str, default = "networks/")
    parser.add_argument("-sgn", "--SingleGeneExpansionList", nargs=1,
                        help="csv file containing gene list downloaded from OpenTarget", 
                        type=str, default = "targets_associated_with_prostate_carcinoma.csv")
    parser.add_argument("-a", "--annotation", help="annotation file from LBM1819",
                        type=str, default="hgnc_cc_zero_filtered_anno.csv")
    parser.add_argument("-r", "--csvboinc", help="csv file from https://gene.disi.unitn.it/test/gene_h.php",
                        type=str, default="genehome PC-IM history.csv")
    parser.add_argument("-n", "--numbgenes", help="number of transcripts to run on NESSRA", 
                        type=int, default=50)

    # Parses and returns the object containing the params
    return parser.parse_args()


def loader_sgn():

    global gene_list

    with open("./output_folder/already-runned-sgn.csv","r") as file_runned:
        file_runned.readline()
    
        for line in file_runned:
            el = line.split(",")
            gene_list[el[0].strip()] = [el[2],el[3:]]

def loader_lgn():

    global gene_list

    with open("./output_folder/already-runned.csv","r") as file_runned:
        file_runned.readline()
    
        for line in file_runned:
            el = line.split(",")
            gene_list[el[0].strip()] = [el[2],el[3:]]
    

def downl_expansion():

    global gene_list

    count = 0
    for expansion in gene_list:
        count = count + 1
        bashcommand = "rclone copy --drive-shared-with-me remote:experiments_results/" + str(expansion) + "_Hs.expansion remote:"
        print(bashcommand)
        os.system(bashcommand)
        print("{}/{} Done 👍".format(count, len(gene_list)))


if __name__ == "__main__":
    
    gene_list = {}

    args = parse_args()

    os.system("python mapper.py -sgn " + args.SingleGeneExpansionList)    
    loader_sgn()
    os.system("python mapper.py -lgn " + args.LocalGeneNetwork)    
    loader_lgn()

    downl_expansion()
