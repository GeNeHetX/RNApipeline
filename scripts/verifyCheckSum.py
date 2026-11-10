import os
import hashlib
import pandas as pd
from glob import glob
from pathlib import Path
import argparse
from tqdm import tqdm
import warnings

warnings.simplefilter("always", Warning)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project_dir",
        "-dir",
        type=Path,
        required=True,
        help="Path to the project directory containing the fastq.gz files to verify.",
        # default="/home/bioinfo/bob/genomic/01_RawRNASeq/PDAC_FFX_NEOADJUVANT/01_RawFiles"
    )
    parser.add_argument(
        "--md5_path",
        "-md5",
        type=Path,
        required=True,
        help="Path to the md5 csv file containing the md5 hash of the fastq.gz files to compare to.",
        # default="/home/bioinfo/wbobgenomic/01_RawRNASeq/PDAC_FFX_NEOADJUVANT/04_Metadata/md5Tab.tsv"
    )
    return parser.parse_args()


def hash_file(filename, blocksize=2**20):
    m = hashlib.md5()
    with filename.open("rb") as f:
        while True:
            # buf = f.read(blocksize)
            buf = f.read()
            if not buf:
                break
            m.update(buf)

    return m.hexdigest()


def fastq_to_md5(filename):
    if not (filename[-5:] != "fastq" or filename[:-8] != "fastq.gz"):
        raise Exception(f"File must be a fastq file, got : {filename} .")
    return filename.split("/")[-1].split(".")[0] + ".md5"


def get_fastq_files(project_dir):
    print(f"Looking for fastq files in {project_dir}")
    list_file = glob(project_dir)
    list_file = [Path(f).resolve() for f in list_file if "fastq" in f]
    print(f"Found {len(list_file)} files to verify.")
    return list_file


def main(args):
    df_saved = pd.read_csv(args.md5_path, sep="\t")
    print(f"Loaded previous md5 hash file containing {len(df_saved)} samples.")

    list_file = get_fastq_files(os.path.join(args.project_dir, "**"))

    if len(list_file) != len(df_saved):
        # print the files that are not in the md5 file
        not_present = [
            f
            for f in list_file
            if str(f).split("/")[-1] not in df_saved["ID_file"].values
        ]
        print(f"Files not present in the md5 file : \n {not_present} \n \n")
        warnings.warn(
            f"Number of files to process ({len(list_file)}) does not match number of files in md5 file ({len(df_saved)})."
        )

    error_count = 0
    for filename in tqdm(
        list_file, desc="Checking files", unit="files", total=len(list_file)
    ):
        curr_name = str(filename).split("/")[-1]
        if curr_name not in df_saved["ID_file"].values:
            print(f"File {curr_name} is not present in the md5 file.")
            warnings.warn(f"File {curr_name} is not present in the md5 file.")
            error_count += 1
            continue

        hash_key = hash_file(filename)
        true_hash_key = df_saved[df_saved["ID_file"] == curr_name]["Hash"].values[0]
        if hash_key != true_hash_key:
            print(f"Hash key for {curr_name} is not the same as in the md5 file.")
            print(f"Hash key for {curr_name} : {hash_key}")
            print(f"Hash key in md5 file : {true_hash_key}")
            error_count += 1
            raise Warning(
                f"Hash key for {curr_name} is not the same as in the md5 file."
            )

    if error_count > 0:
        print(
            f"Number of correct files : {len(list_file) - error_count}. Number of errors : {error_count}"
        )
        with open("status_md5.txt", "w") as f:
            f.write("ERROR!\n")
            f.write(f"Number of files to process in directory:{len(list_file)}\n")
            f.write(f"Number of files in md5 :{len(df_saved)}\n")
            f.write(f"Number of correct files : {len(list_file) - error_count}. Number of errors : {error_count}\n")
        return 1
    else:
        print("All files have been verified with no errors.")
        with open("status_md5.txt", "w") as f:
            f.write("OK")
        return 0


if __name__ == "__main__":
    args = parse_args()
    main(args)
