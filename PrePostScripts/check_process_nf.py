from pathlib import Path
import pandas as pd
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "-path_res", 
        type=Path,
        required=True,
        help="Path to the result process folder",
    )
    parser.add_argument(
        "-path_csv", 
        type=Path, 
        required=True,
        help="Path to the sample csv file", 
    )
    parser.add_argument(
        "-run_number", 
        type=str, 
        help="Run number")
    return parser.parse_args()


def clean_name(ID_file,suffixes):
    for suff in suffixes:
        if ID_file.endswith(suff):
            return ID_file[:-len(suff)]
    return ID_file

def main(args):
    df_csv = pd.read_csv(args.path_csv,sep=',')
    print(f"Checking {len(df_csv)} samples in {args.path_csv}")

    path_starLog = args.path_res / "StarLog_output"
    process_files = list(path_starLog.glob("**/*StarOutLog.final.out"))
    process_files = pd.DataFrame(process_files, columns=["ID_sample"])
    #process_files.head()
    process_files["ID_Patient"] = process_files.ID_sample.apply(lambda x: x.name.split(".")[0])
    #suff_list=list(df_csv["suffix1"]) # GESTION SUFFIX AU BESOIN
    #process_files.ID_Patient = process_files.ID_Patient.apply(lambda x : clean_name(x,suff_list))
    process_files.ID_Patient = process_files.ID_Patient.apply(lambda x: x.replace("StarOutLog", ""))
    non_processed_samples = df_csv[~df_csv.ID_Sample.isin(process_files.ID_Patient)]

    if len(non_processed_samples) > 0:
        with open("status_process.txt", "w") as f:
            f.write("ERROR!\n")
            f.write(f"WARNING: {len(non_processed_samples)} samples not processed on {len(df_csv)}")
            for sample in non_processed_samples.ID_Sample:
                f.write(f"Sample non processed : {sample}")
        return 1
    else:
        print("All samples processed.")
        with open("status_process.txt", "w") as f:
            f.write("OK")
        return 0


if __name__ == "__main__":
    args = parse_args()
    main(args)
