import os
import pandas as pd
import numpy as np
from sets import Set
import fnmatch
import argparse

script_dirpath = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--data-dirpaths', nargs='+', default=[], required=True, help="Dirpath(s) to output data of MG-CFD runs")
parser.add_argument('--output-dirpath', required=True, help="Dirpath to generated processed data")
args = parser.parse_args()
mg_cfd_output_dirpaths = args.data_dirpaths
prepared_output_dirpath = args.output_dirpath

import imp
imp.load_source('data_utils', os.path.join(script_dirpath, "data-utils.py"))
from data_utils import *

def infer_nranks(op2_csv_filepath):
    d = clean_pd_read_csv(op2_csv_filepath)
    return d["nranks"][0]

def infer_partitioner(slurm_filepath):
    with open(slurm_filepath, "r") as f:
        for line in f:
            if line.startswith("partitioner"):
                partitioner = line.replace('\n','').split('=')[1]
                return partitioner
    return ""

def collate_csvs():
    cats = ["PerfData", "PAPI", "op2_performance_data"]

    dirpaths = mg_cfd_output_dirpaths

    for cat in cats:
        print("Collating " + cat)
        df_all = None

        for mg_cfd_output_dirpath in dirpaths:
            sub_dirnames = [i for i in os.listdir(mg_cfd_output_dirpath) if os.path.isdir(os.path.join(mg_cfd_output_dirpath, i))]
            if len(sub_dirnames) == 0:
                raise Exception("No subfolders detected in job directory: " + mg_cfd_output_dirpath)

            for d in sub_dirnames:
                dp = os.path.join(mg_cfd_output_dirpath, d)

                for run_root, run_dirnames, run_filenames in os.walk(dp):
                    nranks = -1
                    ## Need to infer partitioner for op2_performance_data.csv:
                    partitioner = ""
                    batch_fp = ""
                    for f in run_filenames:
                        if f.endswith(".batch") or (f.startswith("run") and f.endswith(".sh")):
                            batch_fp = os.path.join(dp, f)
                            partitioner = infer_partitioner(batch_fp)

                    if batch_fp == "":
                        raise Exception("Failed to find batch file in: " + dp)
                    if partitioner == "":
                        raise Exception("Failed to infer partitioner in: " + batch_fp)
                    op2_perf_data_filepath = os.path.join(dp, "op2_performance_data.csv")
                    if not os.path.isfile(op2_perf_data_filepath):
                        print("WARNING: Job {0} failed, skipping".format(d))
                        continue
                    nranks = infer_nranks(op2_perf_data_filepath)
                    if nranks == -1:
                        raise Exception("Failed to infer nranks in: " + op2_perf_data_filepath)

                    for f in run_filenames:
                        if f.endswith("."+cat+".csv") or (f == cat+".csv"):
                            df_filepath = os.path.join(dp, f)
                            df = clean_pd_read_csv(df_filepath)
                            df["nranks"] = nranks
                            if not "partitioner" in df.columns.values:
                                df["partitioner"] = partitioner
                            if df_all is None:
                                df_all = df
                            else:
                                df_all = df_all.append(df, sort=True)
                                # df_all = df_all.append(df)

        if not df_all is None:
            agg_fp = os.path.join(prepared_output_dirpath,cat+".csv")
            if not os.path.isdir(prepared_output_dirpath):
                os.mkdir(prepared_output_dirpath)
            df_all.to_csv(agg_fp, index=False)

def aggregate():
    for cat in ["PerfData", "PAPI", "op2_performance_data"]:
        print("Aggregating " + cat)
        df_filepath = os.path.join(prepared_output_dirpath,cat+".csv")
        if not os.path.isfile(df_filepath):
            continue
        df = clean_pd_read_csv(df_filepath)
        job_id_colnames = get_job_id_colnames(df)
        data_colnames = get_data_colnames(df)

        ## Compute per-rank average across repeat runs:
        df_agg = df.groupby(job_id_colnames, as_index=False)
        df_mean = df_agg.mean()
        df_mean_out_filepath = os.path.join(prepared_output_dirpath, cat+".mean.csv")
        df_mean.to_csv(df_mean_out_filepath, index=False)

        df_agg = df.groupby(job_id_colnames)
        df_sd = df_agg.std(ddof=0)
        df_sd = df_sd.reset_index()
        df_sd = df_sd.rename(index=str, columns={s:s+"_sd" for s in data_colnames})
        
        df_sd_pct_out_filepath = os.path.join(prepared_output_dirpath, cat+".sd_pct.csv")
        df_sd_pct = df_sd.merge(df_mean)
        for dc in data_colnames:
            df_sd_pct[dc+"_sd_pct"] = 0
            mask = df_sd_pct[dc] != 0
            df_sd_pct.loc[mask, dc+"_sd_pct"] = df_sd_pct.loc[mask, dc+"_sd"] / df_sd_pct.loc[mask, dc]
            df_sd_pct = df_sd_pct.drop(columns=[dc, dc+"_sd"])
        df_sd_pct = df_sd_pct.round(4)
        df_sd_pct.to_csv(df_sd_pct_out_filepath, index=False)

collate_csvs()
aggregate()

print("Aggregated data written to folder: " + prepared_output_dirpath)
