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

flux_data = op2_data.loc[op2_data["loop"]=="flux"]
wall_data = op2_data.loc[op2_data["loop"]=="WALLTIME"]

if "nranks" in op2_data.columns:
	nranks_counts = Set(op2_data["nranks"])

	# # rank vs flux sync:
	# for n in nranks_counts:
	# 	n_data = flux_data[flux_data["nranks"]==n]
	# 	fig = plt.figure(figsize=fig_dims)
	# 	x = "rank"
	# 	y = "sync.percall"
	# 	fig.suptitle("{0} vs {1}".format(x, y))
	# 	ax = fig.add_subplot(1,1,1)
	# 	ax.set_xlabel(x)
	# 	ax.set_ylabel(y)
	# 	loops_ranked_by_sync_cost = n_data[n_data["rank"]==0].sort_values(y, ascending=False)["loop"].values
	# 	for loop in loops_ranked_by_sync_cost:
	# 		loop_times = n_data[n_data["loop"]==loop]
	# 		ax.plot(loop_times[x], loop_times[y], label=loop)
	# 	ax.legend()
	# 	png_filename = 'rank-vs-flux-sync.N={:05d}.png'.format(n)
	# 	png_filepath = os.path.join(output_dirpath, png_filename)
	# 	plt.savefig(png_filepath)
	# 	plt.close(fig)

	# ranks vs flux() sync (mean and stdev):
	flux_data_grp = flux_data.groupby(["nranks"], as_index=False)
	fig = plt.figure(figsize=fig_dims)
	x = "nranks"
	y = "sync.percall"
	fig.suptitle("{0} vs {1}".format(x, y))
	ax = fig.add_subplot(1,1,1)
	ax.set_xlabel(x)
	ax.set_ylabel(y)
	plt.errorbar(flux_data_grp.mean()[x], flux_data_grp.mean()[y], flux_data_grp.std()[y])
	png_filename = 'nranks-vs-flux-sync.png'
	png_filepath = os.path.join(output_dirpath, png_filename)
	plt.savefig(png_filepath)
	plt.close(fig)

	# nranks vs walltime:
	wall_data_grp = wall_data.groupby(["nranks"], as_index=False)
	fig = plt.figure(figsize=fig_dims)
	x = "nranks"
	y = "total time"
	fig.suptitle("{0} vs {1}".format(x, y))
	ax = fig.add_subplot(1,1,1)
	ax.set_xlabel(x)
	ax.set_ylabel(y)
	plt.plot(wall_data_grp.max()[x], wall_data_grp.max()[y])
	ax.set_ylim([0.0, 45.0])
	png_filename = 'nranks-vs-walltime.png'
	png_filepath = os.path.join(output_dirpath, png_filename)
	plt.savefig(png_filepath)
	plt.close(fig)
