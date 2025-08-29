#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np
from pathlib import Path
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project_dir",
        "-p",
        type=str,
        required=True,
        help="Path to the project directory containing the fastq.gz files.",
        # default="/shared/projects/pdacrna/DUMP/"
    )
    parser.add_argument(
        "--sample_csv",
        "-csv",
        type=Path,
        required=True,
        help="Path to the sample csv files to process in the project.",
        # default="/shared/projects/pdacrna/DUMP/project.csv"
    )
    return parser.parse_args()

def extract_suffix(file_name):
    # Trouver la position du dernier '_R'
    last_r_pos = file_name.rfind('_R')
    if last_r_pos == -1:
        return None, None  # Si '_R' n'est pas trouvé dans le nom
    # Extraire tout après le dernier '_R' et avant '.fastq.gz', y compris le '_'
    suffix = file_name[last_r_pos : file_name.rfind('.fastq.gz')]
    ID = file_name[:last_r_pos]
    return ID, suffix

def main(args):
    print("Checking the number of samples in the project directory...")
    df_csv = pd.read_csv(args.sample_csv,sep=',')
    fileList = [f for f in Path(args.project_dir).iterdir() if f.is_file()]
    nameList = [f.name for f in fileList]
    
    # get the suffixes
    dico_suffix = {}
    for file_name in nameList:
        ID, suffix = extract_suffix(file_name)
        if ID is not None and suffix is not None:
            if ID not in dico_suffix:
                dico_suffix[ID] = []
            dico_suffix[ID].append(suffix)
            if len(dico_suffix[ID]) > 1:
                dico_suffix[ID] = sorted(dico_suffix[ID], key=lambda x: (x[:-1], int(x[-1])) if x[-1].isdigit() else x)
        else:
            print(f"File name '{file_name}' does not contain '_R' or '.fastq.gz'.") 
    
    #check if the sample in the csv are in the fileList
    missing_samples = df_csv[~df_csv.ID_Sample.isin(dico_suffix.keys())]
    if len(missing_samples) == 0:
        suffix_df = pd.DataFrame([{'ID_Sample': k, 'suffix1': v[0], 'suffix2': v[1]} for k, v in dico_suffix.items() if len(v) >= 2])
        df_updated = df_csv.merge(suffix_df, on="ID_Sample", how="left")
        df_updated.to_csv("sampleChecked.csv", sep=',', index=False)
        print("All sample found")
        with open("status_samples.txt", "w") as f:
            f.write("OK")
        return 0

    with open("status_samples.txt", "w") as f:
        f.write("ERROR!\n")
        f.write("Check the suffixes or file names for possible errors.\n")
        for index, row in missing_samples.iterrows():
            f.write("Missing sample: " + row['ID_Sample'])
    
    df_csv.to_csv("sampleChecked.csv", sep=',', index=False) # to avoid nextflow to crash
    return 1


if __name__ == "__main__":
    args = parse_args()
    main(args)