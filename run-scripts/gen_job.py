import os, shutil, stat, sys
import json, argparse, re
import itertools
import math
from tempfile import NamedTemporaryFile

script_dirpath = os.path.dirname(os.path.realpath(__file__))
template_dirpath = os.path.join(os.path.dirname(script_dirpath), "run-templates")
app_dirpath = os.path.dirname(script_dirpath)

js_to_filename = {}
js_to_filename[""] = ""
js_to_filename["slurm"] = "slurm.sh"
js_to_filename["moab"] = "moab.sh"
js_to_filename["lsf"] = "lsf.sh"
js_to_filename["pbs"] = "pbs.sh"

js_to_submit_cmd = {}
js_to_submit_cmd[""] = ""
js_to_submit_cmd["slurm"] = "sbatch"
js_to_submit_cmd["moab"] = "msub"
js_to_submit_cmd["lsf"] = "bsub"
js_to_submit_cmd["pbs"] = "qsub"

defaults = {}
# Compilation:
defaults["compiler"] = "intel"
defaults["cpp wrapper"] = ""
defaults["mpicpp wrapper"] = ""
defaults["openmp"] = False
defaults["mpi"] = False
defaults["cuda"] = False
defaults["openacc"] = False
defaults["openmp4"] = False
# Job scheduling:
defaults["unit walltime"] = 0.0
defaults["project code"] = ""
defaults["partition"] = ""
defaults["job scheduler"] = ""
# MG-CFD execution:
defaults["num threads"] = [1]
defaults["num repeats"] = 1
defaults["mg cycles"] = 50

def get_key_value(profile, cat, key):
    if not cat in profile.keys():
        raise Exception("Cat '{0}' not in json.".format(cat))

    if key in profile[cat]:
        return profile[cat][key]
    else:
        if key in defaults:
            return defaults[key]
        else:
            raise Exception("Mandatory key '{0}' not present in cat '{1}' of json".format(key, cat))

def py_sed(filepath, from_str, to_str, delete_line_if_none=False):
    with open(filepath, "r") as f:
        lines = f.readlines()
    with open(filepath, "w") as f:
        for line in lines:
            if from_str in line and to_str == "" and delete_line_if_none:
                # Delete line
                pass
            elif isinstance(to_str, str):
                f.write(re.sub(from_str,     to_str,  line))
            else:
                f.write(re.sub(from_str, str(to_str), line))

def delete_folder_contents(dirpath):
    print("Deleting contents of folder: " + dirpath)
    for f in os.listdir(dirpath):
        fp = os.path.join(dirpath, f)
        if os.path.isdir(fp):
            shutil.rmtree(fp)
        else:
            os.remove(fp)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', required=True)
    args = parser.parse_args()
    with open(args.json, 'r') as f:
        profile = json.load(f)

    ## Read file/folder paths and prepare folders:
    jobs_dir = profile["setup"]["jobs dir"]
    if jobs_dir[0] != '/':
        jobs_dir = os.path.join(app_dirpath, jobs_dir)
    if not os.path.isdir(jobs_dir):
        os.mkdir(jobs_dir)
    else:
        delete_folder_contents(jobs_dir)

    data_dirpath = profile["run"]["data dirpath"]
    if data_dirpath[0] != '/':
        data_dirpath = os.path.join(app_dirpath, data_dirpath)

    ## Read parameters from json:
    job_queue = get_key_value(profile, "setup", "partition")
    project_code = get_key_value(profile, "setup", "project code")
    js = get_key_value(profile, "setup", "job scheduler")

    compiler = get_key_value(profile, "compile", "compiler")
    cpp_wrapper = get_key_value(profile, "compile", "cpp wrapper")
    mpicpp_wrapper = get_key_value(profile, "compile", "mpicpp wrapper")
    use_mpi = get_key_value(profile, "compile", "mpi")
    use_cuda = get_key_value(profile, "compile", "cuda")
    use_openmp = get_key_value(profile, "compile", "openmp")
    use_openacc = get_key_value(profile, "compile", "openacc")
    use_openmp4 = get_key_value(profile, "compile", "openmp4")
    if use_cuda:
        if use_openmp:
            raise Exception("Cannot combine CUDA and OpenMP")
        if use_openmp4:
            raise Exception("Cannot combine CUDA and OpenMP 4")
        if use_openacc:
            raise Exception("Cannot combine CUDA and OpenACC")
    if use_openmp:
        if use_openmp4:
            raise Exception("Cannot combine OpenMP and OpenMP 4")
        if use_openacc:
            raise Exception("Cannot combine OpenMP and OpenACC")
    if use_openmp4:
        if use_mpi:
            raise Exception("Cannot combine OpenMP4 and MPI")
        if use_openacc:
            raise Exception("Cannot combine OpenMP4 and OpenACC")
    if use_openacc:
        if use_mpi:
            raise Exception("Cannot combine OpenACC and MPI")

    if use_mpi:
        if use_cuda:
            bin_filename = "mgcfd_mpi_cuda"
        elif use_openmp:
            bin_filename = "mgcfd_mpi_cuda"
        else:
            bin_filename = "mgcfd_mpi_cuda"

    else:
        if use_cuda:
            bin_filename = "mgcfd_mpi_cuda"
        elif use_openmp:
            bin_filename = "mgcfd_mpi_cuda"
        elif use_openmp4:
            bin_filename = "mgcfd_mpi_cuda"
        elif use_openacc:
            bin_filename = "mgcfd_mpi_cuda"
        else:
            bin_filename = "mgcfd_seq"
    bin_filepath = os.path.join(app_dirpath, "bin", "bin_filename")

    num_repeats = get_key_value(profile, "run", "num repeats")
    mg_cycles = get_key_value(profile, "run", "mg cycles")
    validate_solution = get_key_value(profile, "run", "validate solution")
    mgcfd_unit_runtime_secs = get_key_value(profile, "run", "unit walltime")
    partitioners = get_key_value(profile, "run", "partitioners")
    
    num_nodes_range = get_key_value(profile, "run", "num nodes")
    num_tpn_range = get_key_value(profile, "run", "num tasks per node")
    num_threads_range = get_key_value(profile, "run", "num threads per task")

    num_jobs = len(num_nodes_range) * len(num_tpn_range) * len(num_threads_range) * num_repeats * len(partitioners)

    submit_all_filepath = os.path.join(jobs_dir, "submit_all.sh")
    submit_all_file = open(submit_all_filepath, "w")
    submit_all_file.write("#!/bin/bash\n")
    submit_all_file.write("set -e\n")
    submit_all_file.write("set -u\n")
    submit_all_file.write("\n")
    js_filename = js_to_filename[js]
    submit_all_file.write("# {0}:\n".format(js))
    submit_all_file.write("submit_cmd={0}\n\n".format(js_to_submit_cmd[js]))
    submit_all_file.write("num_jobs={0}\n\n".format(num_jobs))

    with open(os.path.join(jobs_dir, "papi.conf"), "w") as f:
        f.write("PAPI_TOT_INS\n")
        f.write("PAPI_TOT_CYC\n")

    n = 0
    for num_nodes in num_nodes_range:
        for num_tpn in num_tpn_range:
            nt = num_nodes * num_tpn
            for num_thr in num_threads_range:
                for repeat in range(num_repeats):
                    for partitioner in partitioners:
                        n += 1
                        job_id = str(n).zfill(3)

                        print("Creating job {0}/{1}".format(n, num_jobs))

                        job_dir = os.path.join(jobs_dir, job_id)
                        if not os.path.isdir(job_dir):
                            os.mkdir(job_dir)

                        ## Link to papi config file:
                        papi_dest_filepath = os.path.join(job_dir, "papi.conf")
                        if os.path.isfile(papi_dest_filepath):
                            os.remove(papi_dest_filepath)
                        os.symlink(os.path.join(jobs_dir, "papi.conf"), papi_dest_filepath)

                        ## Instantiate MG-CFD run script:
                        job_run_filepath = os.path.join(job_dir, "run-mgcfd.sh")
                        shutil.copyfile(os.path.join(template_dirpath, "run-mgcfd.sh"), job_run_filepath)

                        ## Instantiate job scheduling header:
                        if js != "":
                            js_filepath = os.path.join(job_dir, js_filename)
                            shutil.copyfile(os.path.join(template_dirpath, js_filename), js_filepath)

                        ## Combine into a single batch submission script:
                        if js == "":
                            batch_filename = "run.sh"
                        else:
                            batch_filename = js+".batch"
                        batch_filepath = os.path.join(job_dir, batch_filename)
                        with open(batch_filepath, "w") as f_out:
                            if js != "":
                                with open(js_filepath, "r") as f_in:
                                    for line in f_in.readlines():
                                        f_out.write(line)
                                os.remove(js_filepath)
                            else:
                                f_out.write("#!/bin/bash\n")
                            f_out.write("\n\n")
                            with open(job_run_filepath, "r") as f_in:
                                for line in f_in.readlines():
                                    f_out.write(line)
                            os.remove(job_run_filepath)

                        ## Now replace variables in script:

                        ## - File/dir paths:
                        py_sed(batch_filepath, "<RUN_OUTDIR>", job_dir)
                        py_sed(batch_filepath, "<APP_DIRPATH>", app_dirpath)
                        py_sed(batch_filepath, "<DATA_DIRPATH>", data_dirpath)

                        ## - Scheduling:
                        py_sed(batch_filepath, "<RUN ID>", job_id)
                        py_sed(batch_filepath, "<PARTITION>", job_queue, True)
                        py_sed(batch_filepath, "<PROJECT CODE>", project_code, True)

                        ## - Parallelism:
                        py_sed(batch_filepath, "<TPN>", num_tpn)
                        py_sed(batch_filepath, "<NODES>", num_nodes)
                        py_sed(batch_filepath, "<NTHREADS>", num_thr)
                        py_sed(batch_filepath, "<NTASKS>", nt)
                        py_sed(batch_filepath, "<NCPUS_PER_NODE>", num_tpn*num_thr)
                        py_sed(batch_filepath, "<NTHREADS>", num_thr)

                        ## - Compilation:
                        py_sed(batch_filepath, "<COMPILER>", compiler)
                        py_sed(batch_filepath, "<CPP_WRAPPER>", cpp_wrapper)
                        py_sed(batch_filepath, "<MPICPP_WRAPPER>", mpicpp_wrapper)
                        py_sed(batch_filepath, "<MPI>", str(use_mpi).lower())
                        py_sed(batch_filepath, "<CUDA>", str(use_cuda).lower())
                        py_sed(batch_filepath, "<OPENMP>", str(use_openmp).lower())
                        py_sed(batch_filepath, "<OPENMP4>", str(use_openmp4).lower())
                        py_sed(batch_filepath, "<OPENACC>", str(use_openacc).lower())

                        ## - Execution:
                        py_sed(batch_filepath, "<PARTITIONER>", partitioner)
                        py_sed(batch_filepath, "<MG_CYCLES>", mg_cycles)
                        if validate_solution:
                            py_sed(batch_filepath, "<VALIDATE_SOLUTION>", "true")
                        else:
                            py_sed(batch_filepath, "<VALIDATE_SOLUTION>", "false")

                        ## - Walltime estimation:
                        if mgcfd_unit_runtime_secs == 0.0:
                            est_runtime_hours = 0
                            est_runtime_minutes = 30
                        else:
                            est_runtime_secs = float(mgcfd_unit_runtime_secs*mg_cycles) / math.sqrt(float(nt*num_thr))
                            est_runtime_secs *= 1.2 ## Allow for estimation error
                            est_runtime_secs += 20  ## Add time for file load
                            est_runtime_secs = int(round(est_runtime_secs))
                            est_runtime_hours = est_runtime_secs/60/60
                            est_runtime_secs -= est_runtime_hours*60*60
                            est_runtime_minutes = est_runtime_secs/60
                            est_runtime_secs -= est_runtime_minutes*60
                            if est_runtime_secs > 0:
                                est_runtime_minutes += 1
                                est_runtime_secs = 0
                        py_sed(batch_filepath, "<HOURS>", str(est_runtime_hours).zfill(2))
                        py_sed(batch_filepath, "<MINUTES>", str(est_runtime_minutes).zfill(2))
                        py_sed(batch_filepath, "<RUN_ID>", job_id)

                        ## Make batch script executable:
                        os.chmod(batch_filepath, 0755)

                        ## Add an entry to submit_all.sh:
                        # submit_all_file.write("\n\n")
                        # submit_all_file.write('if [ ! -f "{0}"/Times.csv ]; then\n'.format(job_dir))
                        # submit_all_file.write("  basedir=`pwd`\n")
                        # submit_all_file.write("  cd {0}\n".format(job_dir))
                        # if js == "":
                        #     ## Without a job scheduler to log STDOUT, need to do this manually:
                        #     submit_all_file.write('  echo "Executing job {0}/{1}"\n'.format(n, num_jobs))
                        #     submit_all_file.write('  eval "$submit_cmd" ./{0} > submit.log\n'.format(batch_filename))
                        # else:
                        #     submit_all_file.write('  eval "$submit_cmd" ./{0}\n'.format(batch_filename))
                        # submit_all_file.write('  cd "$basedir"\n')
                        # submit_all_file.write('fi')

                        # Copy template 'submit.sh' to a temp file:
                        submit_tmp = NamedTemporaryFile(prefix='myprefix')
                        if not os.path.isfile(submit_tmp.name):
                            print("NamedTemporaryFile() failed to actually create the file")
                            sys.exit(-1)
                        submit_tmp_filepath = submit_tmp.name
                        with open(os.path.join(template_dirpath, "submit.sh"), 'r+b') as f:
                            shutil.copyfileobj(f, submit_tmp)
                        submit_tmp.seek(0)
                        # Instantiate and append to submit_all.sh: 
                        py_sed(submit_tmp_filepath, "<RUN_DIRPATH>",    job_dir)
                        py_sed(submit_tmp_filepath, "<BATCH_FILENAME>", batch_filename)
                        py_sed(submit_tmp_filepath, "<BIN_FILEPATH>",   bin_filepath)
                        py_sed(submit_tmp_filepath, "<JOB_NUM>", n)
                        py_sed(submit_tmp_filepath, "<NUM_JOBS>", num_jobs)
                        submit_all_file.write("\n\n")
                        with open(submit_tmp_filepath, 'r') as f:
                            for line in f:
                                # if not re.match(r''+"^[\s]*#", line):
                                submit_all_file.write(line)

                        ## Now close (and delete) submit_tmp file
                        submit_tmp.close()
                        if os.path.isfile(submit_tmp_filepath):
                            os.unlink(submit_tmp_filepath)


    submit_all_file.close()
    st = os.stat(submit_all_filepath)
    os.chmod(submit_all_filepath, st.st_mode | stat.S_IEXEC)
