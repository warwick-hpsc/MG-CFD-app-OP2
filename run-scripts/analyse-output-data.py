import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
fig_dims = (10,8)

import os
import argparse
from copy import deepcopy

script_dirpath = os.path.dirname(os.path.realpath(__file__))

import imp
imp.load_source('data_utils', os.path.join(script_dirpath, "data-utils.py"))
from data_utils import *

parser = argparse.ArgumentParser()
parser.add_argument('--data-dirpath', required=True, help="Dirpath to processed data")
parser.add_argument('--output-dirpath', required=False, help="Dirpath to generated processed data")
args = parser.parse_args()
data_dirpath = args.data_dirpath

if not args.output_dirpath is None:
    output_dirpath = args.output_dirpath
else:
    output_dirpath = data_dirpath +".analysed"
if not os.path.isdir(output_dirpath):
	os.mkdir(output_dirpath)

op2_data = pd.read_csv(os.path.join(data_dirpath, "op2_performance_data.mean.csv"), keep_default_na=False)

## Clean dataset:
# Remove non-varying job-id columns:
job_id_colnames = get_job_id_colnames(op2_data)
for j in job_id_colnames:
	vals = Set(op2_data[j])
	if len(vals) == 1:
		op2_data = op2_data.drop(columns=[j])
# Remove unwanted data columns:
unwanted_columns = ["plan time", "GB used", "GB total"]
for col in unwanted_columns:
	if col in op2_data.columns:
		op2_data = op2_data.drop(columns=[col])

job_id_colnames = get_job_id_colnames(op2_data)

## Sum across kernels for overall walltime:
group_col_names = deepcopy(job_id_colnames)
group_col_names.remove("kernel name")
df_agg = op2_data.groupby(group_col_names, as_index=False)
op2_walltimes = df_agg.sum()
op2_walltimes["kernel name"] = "WALLTIME"
op2_data = pd.concat([op2_data, op2_walltimes], sort=False)

# Drop unwanted loops:
wanted_loops = ["compute_flux_edge_kernel", "WALLTIME"]
f = op2_data["kernel name"]=="gibberish-qwerty"
for l in wanted_loops:
	f = np.logical_or(f, op2_data["kernel name"]==l)
op2_data = op2_data[f]
op2_data.loc[:,'kernel name'].replace("compute_flux_edge_kernel", "flux", inplace=True)

# Rename columns:
op2_data = op2_data.rename(index=str, columns={"kernel name":"loop", "mpi time":"sync"})

# Calculate per-call times and counts:
timers = ["total time", "sync"]
for t in timers:
	op2_data[t+".percall"] = op2_data[t] / op2_data["count"]
op2_data = op2_data.drop(columns=["count"])

op2_data.to_csv(os.path.join(output_dirpath, "op2_data.csv"), index=False)

# ## Calculate MPI % (per rank):
# op2_data["sync pct"] = op2_data["sync"] / op2_data["total time"]
# # ... and write out:
# mpi_pct_filepath = os.path.join(output_dirpath, "mpi_pct_rankwise.csv")
# if not os.path.isfile(mpi_pct_filepath):
# 	mpi_pct = op2_data.copy()
# 	unwanted_columns = ["sync", "GB used"]
# 	for col in unwanted_columns:
# 		if col in mpi_pct.columns:
# 			mpi_pct = mpi_pct.drop(columns=[col])
# 	if not os.path.isdir(output_dirpath):
# 		os.mkdir(output_dirpath)
# 	mpi_pct["sync pct"] = np.round(mpi_pct["sync pct"], decimals=4)
# 	mpi_pct.to_csv(mpi_pct_filepath, index=False)

#####################

## Plot:

nranks_counts = Set(op2_data["nranks"])
for n in nranks_counts:
	n_data = op2_data[op2_data["nranks"]==n]
	# rank vs sync:
	fig = plt.figure(figsize=fig_dims)
	x = "rank"
	y = "sync.percall"
	fig.suptitle("{0} vs {1}".format(x, y))
	ax = fig.add_subplot(1,1,1)
	ax.set_xlabel(x)
	ax.set_ylabel(y)
	loops_ranked_by_sync_cost = n_data[n_data["rank"]==0].sort_values(y, ascending=False)["loop"].values
	for loop in loops_ranked_by_sync_cost:
		loop_times = n_data[n_data["loop"]==loop]
		ax.plot(loop_times[x], loop_times[y], label=loop)
	ax.legend()
	# plt.savefig('rank-vs-sync.png')
	png_filename = 'rank-vs-sync.N={:05d}.png'.format(n)
	png_filepath = os.path.join(output_dirpath, png_filename)
	plt.savefig(png_filepath)
	plt.close(fig)

#####################

# ## Now perform statistics across ranks of each run (worst case, median, sum etc):

# group_col_names = deepcopy(job_id_colnames)
# group_col_names.remove("rank")

# op2_data_fluxes = op2_data[op2_data["kernel name"]=="compute_flux_edge_kernel"]
# op2_data_wall = op2_data[op2_data["kernel name"]=="WALLTIME"]

# op2_data_fluxes_grp = op2_data_fluxes.groupby(group_col_names, as_index=False)
# op2_data_wall_grp = op2_data_wall.groupby(group_col_names, as_index=False)

# # Worst case:
# op2_data_mpi_wc_filepath = os.path.join(output_dirpath, "mpi_wc_runwise.csv")
# if not os.path.isfile(op2_data_mpi_wc_filepath):
# 	op2_data_fluxes_mpi_wc = op2_data_fluxes.loc[op2_data_fluxes_grp["sync"].idxmax()]
# 	op2_data_wall_mpi_wc = op2_data_wall.loc[op2_data_wall_grp["sync"].idxmax()]
# 	op2_data_mpi_wc = pd.concat([op2_data_fluxes_mpi_wc, op2_data_wall_mpi_wc])
# 	op2_data_mpi_wc.to_csv(op2_data_mpi_wc_filepath, index=False)

# # Sum across ranks:
# op2_data_fluxes_mpi_sum = op2_data_fluxes_grp.sum()
# op2_data_wall_mpi_sum = op2_data_wall_grp.sum()
# op2_data_mpi_sum = pd.concat([op2_data_fluxes_mpi_sum, op2_data_wall_mpi_sum])
# op2_data_mpi_sum["sync pct"] = op2_data_mpi_sum["sync"] / op2_data_mpi_sum["total time"]
# op2_data_mpi_sum["sync pct"] = np.round(op2_data_mpi_sum["sync pct"], decimals=4)
# op2_data_mpi_sum_filepath = os.path.join(output_dirpath, "mpi_sum_runwise.csv")
# if not os.path.isfile(op2_data_mpi_sum_filepath):
# 	unwanted_columns = ["GB used"]
# 	for col in unwanted_columns:
# 		if col in op2_data_mpi_sum.columns:
# 			op2_data_mpi_sum = op2_data_mpi_sum.drop(columns=[col])
# 	op2_data_mpi_sum.to_csv(op2_data_mpi_sum_filepath, index=False)

# # Median across ranks:
# op2_data_fluxes_mpi_median = op2_data_fluxes_grp.median()
# op2_data_wall_mpi_median = op2_data_wall_grp.median()
# op2_data_mpi_median = pd.concat([op2_data_fluxes_mpi_median, op2_data_wall_mpi_median])
# op2_data_mpi_median["sync pct"] = op2_data_mpi_median["sync"] / op2_data_mpi_median["total time"]
# op2_data_mpi_median["sync pct"] = np.round(op2_data_mpi_median["sync pct"], decimals=4)
# op2_data_mpi_median_filepath = os.path.join(output_dirpath, "mpi_median_runwise.csv")
# if not os.path.isfile(op2_data_mpi_median_filepath):
# 	unwanted_columns = ["GB used"]
# 	for col in unwanted_columns:
# 		if col in op2_data_mpi_median.columns:
# 			op2_data_mpi_median = op2_data_mpi_median.drop(columns=[col])
# 	op2_data_mpi_median.to_csv(op2_data_mpi_median_filepath, index=False)
