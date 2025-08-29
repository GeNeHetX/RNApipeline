import os
import hashlib
import pandas as pd
from glob import glob
from pathlib import Path
import argparse
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project_dir",
        "-p",
        type=Path,
        required=True,
        help="Path to the project directory containing the files to hash.",
    )
    parser.add_argument(
        "--sample_csv",
        "-c",
        type=Path,
        required=True,
        help="Path to csv file containing the sample ID and suffixes.",
    )
    parser.add_argument(
        "--md5_path",
        "-md5",
        type=Path,
        required=True,
        help="Path to the md5 tsv file containing the md5 hash of the fastq.gz files. If the file does not exist, it will be created.",
        default="/home/bioinfo/wbobgenomic/01_RawRNASeq/PDAC_FFX_NEOADJUVANT/03_QC/md5Tab.tsv",
    )
    parser.add_argument(
        "--extension",
        "-e",
        type=str,
        required=False,
        nargs="+",
        help="Extension of the files to process. Can be multiple.",
        default=["fastq.gz", "fastq"],
    )
    return parser.parse_args()


def hash_file(filename, blocksize=2**20):
    
    m = hashlib.md5()
    with filename.open("rb") as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update(buf)

    return m.hexdigest()


def filename_in_extension(filename, extension):
    return any([str(filename).endswith(ext) for ext in extension])


def get_files(project_dir, df_samples, extensions):
    print(f"Looking for files in {project_dir}")
    # Générer les noms attendus
    expected_filenames = set()
    for i, row in df_samples.iterrows():
        for ext in extensions:
            expected_filenames.add(f"{row['ID_Sample']}{row['suffix1']}.{ext}")
            expected_filenames.add(f"{row['ID_Sample']}{row['suffix2']}.{ext}")
    list_file = []
    for file in project_dir.glob(f"**/*"):
        if file.name in expected_filenames:
            list_file.append(file)
    print(f"Found {len(list_file)} files to process.")
    return list_file


def main(args):
    
    try:
        df_saved = pd.read_csv(args.md5_path, sep="\t")
        print(f"Loaded previous md5 hash file containing {len(df_saved)} samples.")
    except FileNotFoundError:
        df_saved = pd.DataFrame(columns=["ID_file", "Hash", "Path"])
        print(f"File {args.md5_path} not found, creating new md5 file.")

    # df_saved = pd.DataFrame(columns=["ID_file", "Hash", "Path"])
    df_samples = pd.read_csv(args.sample_csv, sep=",")
    list_file = get_files(args.project_dir, df_samples, args.extension)
    # Keep only new samples
    list_file = [f for f in list_file if f.name not in df_saved["ID_file"].values]
    print(f"Found {len(list_file)} new files to process.")

    dict_list = []
    for filename in tqdm(list_file, desc="Hashing files", unit="files", total=len(list_file)):
        curr_name = filename.name
        # get hash key
        hash_key = hash_file(filename)
        new_dict = {"ID_file": curr_name, "Hash": hash_key, "Path": filename}
        dict_list.append(new_dict)

        new_df = pd.concat([df_saved, pd.DataFrame(dict_list)])
        new_df.sort_values(by=["ID_file"], inplace=True)
        new_df.to_csv(args.md5_path, sep="\t", index=False)
        # Save at each step in case of crash

    print("Done !")
    with open("status_md5.txt", "w") as f:
        f.write("OK")


if __name__ == "__main__":
    args = parse_args()
    main(args)
