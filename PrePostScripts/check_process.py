from pathlib import Path
import pandas as pd
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "-path_process", 
        type=Path,
        help="Path to the process folder", 
        default="/shared/projects/pdacrna/DUMP",
    )
    parser.add_argument(
        "-path_md5", 
        type=Path, 
        help="Path to the md5 file", 
        default="/shared/projects/pdacrna/DUMP/md5Tab.tsv",
    )
    parser.add_argument(
        "--suffixes",
        "-s",
        nargs='+',
        type=str)
    
    parser.add_argument(
        "-run_number", 
        type=str, 
        help="Run number")
    # En vrai je suis pas sur du tout, peut-être tu peux avoir une erreur dans le pipeline sur le R2 
    # sans que le fichier de log ne soit impacté...
    return parser.parse_args()


def clean_name(ID_file,suffixes):
    for suff in suffixes:
        if ID_file.endswith(suff):
            return ID_file[:-len(suff)]
    return ID_file

def main(args):
    df_md5 = pd.read_csv(args.path_md5, sep="\t")
    df_md5["ID_Patient"] = df_md5["ID_file"].apply(lambda x: x.split(".")[0])
    df_md5["ID_Patient"] = df_md5["ID_Patient"].apply(lambda x : clean_name(x,args.suffixes))

    path_files = args.path_process / args.run_number / "FeatureCounts_output"
    process_files = list(path_files.glob("**/*StarOutLog.final.out"))
    process_files = pd.DataFrame(process_files, columns=["ID_sample"])
    process_files["ID_Patient"] = process_files.ID_sample.apply(lambda x: x.name.split(".")[0])
    process_files.ID_Patient = process_files.ID_Patient.apply(lambda x : clean_name(x,args.suffixes))
    process_files.ID_Patient = process_files.ID_Patient.apply(lambda x: x.replace("StarOutLog", ""))
    non_processed_samples = df_md5[~df_md5.ID_Patient.isin(process_files.ID_Patient)]

    if len(non_processed_samples) > 0:
        print(f"WARNING: {len(non_processed_samples)} samples not processed on {len(df_md5)}")
        for sample in non_processed_samples.ID_Patient:
            print(f"Sample non processed : {sample}")
        print(f"WARNING: {len(non_processed_samples)} samples not processed on {len(df_md5)}")
        with open("status_process.txt", "w") as f:
            f.write("ERROR")
        return 1
    else:
        print("All samples processed.")
        with open("status_process.txt", "w") as f:
            f.write("OK")
        return 0


if __name__ == "__main__":
    args = parse_args()
    main(args)
