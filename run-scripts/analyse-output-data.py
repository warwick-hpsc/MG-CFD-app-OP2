import pandas as pd
import numpy as np
import os
import argparse
from copy import deepcopy

import matplotlib.pyplot as plt
fig_dims = (10,8)

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

## Read in data
op2_data = pd.read_csv(os.path.join(data_dirpath, "op2_performance_data.mean.csv"), keep_default_na=False)
op2_data["loop time"] = op2_data["total time"] - op2_data["mpi time"]
# op2_data.drop(columns=["total time"], inplace=True)

perfdata_filepath = os.path.join(data_dirpath, "PerfData.mean.csv")
if not os.path.isfile(perfdata_filepath):
	flux_iters_df = None
else:
	perfdata_df = pd.read_csv(perfdata_filepath)
	flux_f = perfdata_df["kernel"]=="compute_flux_edge_kernel"
	flux_iters_df = perfdata_df[flux_f][["partitioner", "rank", "nranks", "iters"]]
	flux_iters_df = flux_iters_df.rename(index=str, columns={"iters":"flux_iters"})
	op2_data = op2_data.merge(flux_iters_df, validate="many_to_one")

papidata_filepath = os.path.join(data_dirpath, "PAPI.mean.csv")
if not os.path.isfile(papidata_filepath):
	papi_df = None
else:
	papi_df = pd.read_csv(papidata_filepath)
	flux_f = papi_df["kernel"]=="compute_flux_edge_kernel"
	flux_papi_df = papi_df[flux_f][["PAPI counter", "partitioner", "rank", "thread", "nranks", "count"]]
	offcore_dram_read_events = [
		"OFFCORE_REQUESTS:ALL_DATA_READ",
		"OFFCORE_RESPONSE_[01]:(DMND_DATA_RD|PF_DATA_RD|ANY_DATA):(DMND_DATA_RD|PF_DATA_RD|ANY_DATA)",
		"OFFCORE_RESPONSE_[01]:(DMND_DATA_RD|PF_DATA_RD|ANY_DATA)"
	]
	f = None
	for e in offcore_dram_read_events:
		g = flux_papi_df["PAPI counter"].str.match(r''+e)
		if f is None:
			f = g
		else:
			f = np.logical_or(f, g)
	flux_gb_df = flux_papi_df[f]
	if flux_gb_df.shape[0] > 0:
		flux_gb_df.drop(columns=["PAPI counter"], inplace=True)
		flux_gb_df["GB read"] = flux_gb_df["count"]*64/1e9
		flux_gb_df.drop(columns=["count"], inplace=True)
		op2_data = op2_data.merge(flux_gb_df, validate="many_to_one")

## Clean dataset:
# Remove non-varying job-id columns:
job_id_colnames = get_job_id_colnames(op2_data)
for j in job_id_colnames:
	vals = Set(op2_data[j])
	if len(vals) == 1:
		if j not in ["nranks", "rank"]:
			op2_data = op2_data.drop(columns=[j])
# Remove unwanted data columns:
unwanted_columns = ["plan time", "GB used", "GB total"]
for col in unwanted_columns:
	if col in op2_data.columns:
		op2_data = op2_data.drop(columns=[col])

job_id_colnames = get_job_id_colnames(op2_data)

## Sum across kernels for walltime:
group_col_names = deepcopy(job_id_colnames)
group_col_names.remove("kernel")
df_agg = op2_data.groupby(group_col_names, as_index=False)
op2_walltimes = df_agg.sum()
op2_walltimes["kernel"] = "WALLTIME"
op2_data = pd.concat([op2_data, op2_walltimes], sort=False)

# Drop unwanted loops:
wanted_loops = ["compute_flux_edge_kernel", "unstructured_stream_kernel", "WALLTIME"]
f = op2_data["kernel"].isin(wanted_loops)
op2_data = op2_data[f]
op2_data.loc[:,"kernel"].replace("compute_flux_edge_kernel", "flux", inplace=True)
op2_data.loc[:,"kernel"].replace("unstructured_stream_kernel", "ustream", inplace=True)

# Rename columns:
op2_data = op2_data.rename(index=str, columns={"kernel":"loop", "mpi time":"sync"})

# Calculate per-call times and counts:
timers = ["loop time", "sync"]
for t in timers:
	op2_data[t+".percall"] = (op2_data[t] / op2_data["count"]).round(8)
# op2_data = op2_data.drop(columns=["count"])

bw_df = None
if "flux_iters" in op2_data.columns.values:
	# Calculate mesh bytes transferred:
	# Need to know number of nodes in mesh, not currently exported by MG-CFD or OP2.
	# Rather than code this, can instead infer it as only two meshes are currently 
	# supported. 
	# Determine whether M6-wing or some variant of Rotor37 to assume ratio of edges to nodes, 
	# then apply to 'flux_iters' which is number of edges.
	m6wing_nedges_MG_avg = 609529.25
	op2_data["nedges"] = op2_data["flux_iters"]
	op2_data["nedges/call"] = op2_data["nedges"] / op2_data["count"]
	group_cols = [c for c in ["loop", "nranks"] if c in op2_data.columns.values]
	op2_data_grp = op2_data.groupby(group_cols, as_index=False)
	op2_data_grp_sum = op2_data_grp.sum()[["loop", "nranks", "nedges/call"]]
	nedges = op2_data_grp_sum[op2_data_grp_sum["loop"]=="flux"].reset_index().loc[0,"nedges/call"]
	## In one multigrid V-cycle, the ratio between total numbers of nodes and edges is:
	if abs((nedges-m6wing_nedges_MG_avg)/m6wing_nedges_MG_avg) < 0.1:
		# Assume M6 wing
		print("- note: assuming M6-wing mesh to infer number of nodes")
		edges_per_node = 4.04
	else:
		# Assume Rotor37 mesh
		print("- note: assuming Rotor37 mesh to infer number of nodes")
		edges_per_node = 3.68

	op2_data["nnodes"] = (op2_data["nedges"]/edges_per_node).round()
	op2_data["edges effective GB read"]  = ((op2_data["nedges"]*3   *8)/1e9).round(5)
	op2_data["nodes effective GB read"]  = ((op2_data["nnodes"]*5*2 *8)/1e9).round(5)
	op2_data["nodes effective GB write"] = ((op2_data["nnodes"]*5   *8)/1e9).round(5)
	op2_data["effective GB read"] = op2_data["edges effective GB read"] + op2_data["nodes effective GB read"]
	op2_data["effective GB write"] = op2_data["nodes effective GB write"]
	op2_data["effective GB"] = op2_data["effective GB read"] + op2_data["effective GB write"]

	group_cols = [c for c in ["loop", "nranks"] if c in op2_data.columns.values]
	op2_data_grp = op2_data.groupby(group_cols, as_index=False)
	op2_data_grp_mean = op2_data_grp.mean()[["loop", "nranks", "loop time"]]
	op2_data_grp_sum = op2_data_grp.sum()[["loop", "nranks", "effective GB", "effective GB read"]]
	gb_sec_data = op2_data_grp_mean.merge(op2_data_grp_sum, validate="one_to_one")
	gb_sec_data["effective read GB/sec"] = (gb_sec_data["effective GB read"] / gb_sec_data["loop time"]).round(1)
	gb_sec_data["effective GB/sec"] = (gb_sec_data["effective GB"] / gb_sec_data["loop time"]).round(1)
	gb_sec_data.drop(columns=["effective GB", "effective GB read", "loop time"], inplace=True)
	bw_df = gb_sec_data

	if "GB read" in op2_data.columns.values:
		# op2_data["actual GB"] = op2_data["actual GB read"]
		## PAPI counter for recording actual DRAM traffic can only capture reads, not writes.
		## But there is a solution.
		##
		## First, note that:
		## actual GB read = effective GB read + wasted GB read
		## wasted GB read = cacheline bytes never read when line is later evicted
		##
		## For edge data, accessed directly with unit stride, there should be zero wasted GB read, 
		## so effective == actual.
		##
		## Then, note that:
		## actual GB read = actual edges GB read + actual nodes GB read
		##
		## This gives us 'actual nodes GB read' without having to calculate ratio of cache misses etc.
		##
		## For each node, 10 doubles are read and 5 written, so can estimate 'actual nodes GB write':
		op2_data["nodes GB read"] = op2_data["GB read"] - op2_data["edges effective GB read"]
		op2_data["nodes estimated GB write"] = op2_data["nodes GB read"] / 2.0
		# op2_data["nodes estimated GB write"] = op2_data["nodes effective GB write"]
		op2_data["estimated GB"] = op2_data["edges effective GB read"] + op2_data["nodes GB read"] + op2_data["nodes estimated GB write"]

		group_cols = [c for c in ["loop", "nranks"] if c in op2_data.columns.values]
		op2_data_grp = op2_data.groupby(group_cols, as_index=False)
		op2_data_grp_mean = op2_data_grp.mean()[["loop", "nranks", "loop time"]]
		op2_data_grp_sum = op2_data_grp.sum()[["loop", "nranks", "GB read", "estimated GB"]]
		gb_sec_data = op2_data_grp_mean.merge(op2_data_grp_sum, validate="one_to_one")
		gb_sec_data["read GB/sec"] = (gb_sec_data["GB read"] / gb_sec_data["loop time"]).round(1)
		gb_sec_data["estimated GB/sec"] = (gb_sec_data["estimated GB"] / gb_sec_data["loop time"]).round(1)
		gb_sec_data.drop(columns=["GB read", "estimated GB", "loop time"], inplace=True)
		bw_df = bw_df.merge(gb_sec_data, validate="one_to_one")

	op2_data.drop(columns=["flux_iters", "nedges", "nnodes", "nedges/call"], inplace=True, errors='ignore')
	op2_data.drop(columns=["edges effective GB read", "nodes effective GB read", "nodes effective GB write"], inplace=True, errors='ignore')
	op2_data.drop(columns=["effective GB read", "effective GB write", "effective GB"], inplace=True, errors='ignore')
	op2_data.drop(columns=["nodes GB read", "nodes estimated GB write"], inplace=True, errors='ignore')
	op2_data.drop(columns=["GB read", "estimated GB"], inplace=True, errors='ignore')

	wanted_loops = ["flux", "ustream"]
	f = bw_df["loop"].isin(wanted_loops)
	bw_df = bw_df[f]

op2_data.to_csv(os.path.join(output_dirpath, "op2_data.csv"), index=False)
if not bw_df is None:
	bw_df.to_csv(os.path.join(output_dirpath, "node_bw_performance.csv"), index=False)

## Calculate MPI % (per rank):
op2_data["sync pct"] = op2_data["sync"] / (op2_data["loop time"] + op2_data["sync"])
# ... and write out:
mpi_pct_filepath = os.path.join(output_dirpath, "mpi_pct_rankwise.csv")
if not os.path.isfile(mpi_pct_filepath):
	mpi_pct = op2_data.copy()
	unwanted_columns = ["sync", "GB used"]
	for col in unwanted_columns:
		if col in mpi_pct.columns:
			mpi_pct = mpi_pct.drop(columns=[col])
	if not os.path.isdir(output_dirpath):
		os.mkdir(output_dirpath)
	mpi_pct["sync pct"] = np.round(mpi_pct["sync pct"], decimals=4)
	mpi_pct.to_csv(mpi_pct_filepath, index=False)

#####################

## Plot:

flux_data = op2_data.loc[op2_data["loop"]=="flux"]
ustream_data = op2_data.loc[op2_data["loop"]=="ustream"]
wall_data = op2_data.loc[op2_data["loop"]=="WALLTIME"].copy()
wall_data["total time"] = wall_data["loop time"] + wall_data["sync"]

cols_to_keep = [c for c in ["nranks", "rank", "thread"] if c in flux_data.columns.values]
flux_with_ustream_data = flux_data.rename(index=str, columns={"loop time":"flux time"})[cols_to_keep+["flux time"]]
flux_with_ustream_data = flux_with_ustream_data.merge(ustream_data.rename(index=str, columns={"loop time":"ustream time"})[cols_to_keep+["ustream time"]])
flux_with_ustream_data["flux_pct_ustream"] = flux_with_ustream_data["flux time"] / flux_with_ustream_data["ustream time"]

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

	if wall_data.shape[0] > 0:
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
		ax.set_ylim([0.0, wall_data[y].max()])
		png_filename = 'nranks-vs-walltime.png'
		png_filepath = os.path.join(output_dirpath, png_filename)
		plt.savefig(png_filepath)
		plt.close(fig)

	if ustream_data.shape[0] > 0:
		# contrast flux() and ustream() scaling:
		flux_data_grp = flux_data.groupby(["nranks"], as_index=False)
		ustream_data_grp = ustream_data.groupby(["nranks"], as_index=False)
		fig = plt.figure(figsize=fig_dims)
		x = "nranks"
		y = "loop time"
		fig.suptitle("{0} vs {1}".format(x, y))
		ax = fig.add_subplot(1,1,1)
		ax.set_xlabel(x)
		ax.set_ylabel(y)
		ax.set_ylim([0.0, flux_data[y].max()])
		ax.set_xlim([flux_data[x].min()-8, flux_data[x].max()+16])
		plt.xticks(np.arange(flux_data[x].min(), flux_data[x].max()+16, 16))
		plt.plot(flux_data_grp.mean()[x], flux_data_grp.mean()[y], 'o-', label="compute_flux()")
		plt.plot(ustream_data_grp.mean()[x], ustream_data_grp.mean()[y], 'o-', label="data bound")
		plt.legend()
		png_filename = 'flux-scaling-vs-mem-bound.png'
		png_filepath = os.path.join(output_dirpath, png_filename)
		plt.savefig(png_filepath)
		plt.close(fig)

		# contrast flux() and ustream() scaling:
		flux_with_ustream_data_grp = flux_with_ustream_data.groupby(["nranks"], as_index=False)
		fig = plt.figure(figsize=fig_dims)
		x = "nranks"
		y = "flux_pct_ustream"
		ylabel = "compute_flux() / unstructured_stream()"
		title = "scaling of " + ylabel
		fig.suptitle(title)
		ax = fig.add_subplot(1,1,1)
		ax.set_xlabel(x)
		ax.set_ylabel(ylabel)
		ax.set_ylim([1.0, 4.0])
		ax.set_xlim([flux_with_ustream_data[x].min()-8, flux_with_ustream_data[x].max()+16])
		plt.xticks(np.arange(flux_with_ustream_data[x].min(), flux_with_ustream_data[x].max()+16, 16))
		plt.plot(flux_with_ustream_data_grp.mean()[x], flux_with_ustream_data_grp.mean()[y], 'o-')
		png_filename = 'flux-scaling-pct-mem-bound.png'
		png_filepath = os.path.join(output_dirpath, png_filename)
		plt.savefig(png_filepath)
		plt.close(fig)
