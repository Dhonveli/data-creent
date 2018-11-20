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

if __name__ == "__main__":
    
    gene_list = {}

    args = parse_args()

    if not os.path.isfile("./output_folder/already-runned.csv"):
        print("You have to run mapper2 before!") 
    else:
        with open("./output_folder/already-runned.csv","r") as file_runned:
            with open("downloaded.csv","r+") as file_down:
                file_runned.readline()
                lines = [line.strip() for line in file_down]   
                count = 0
                for line in file_runned:
                    tmp = line.split(",")
                    if tmp[0] in lines:
                        continue
                    else:
                        count += 1
                        bashcommand = "rclone copy --drive-shared-with-me gDrive:experiments_results/" + str(tmp[0]) + "_Hs.expansion gDrive:"
                        print(bashcommand)
                        os.system(bashcommand)
                        print("{} Done 👍".format(count))
                        file_down.write("%s\n" % tmp[0])
