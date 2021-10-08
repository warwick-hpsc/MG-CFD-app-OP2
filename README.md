MG-CFD OP2
==========================================

OP2 port of [MG-CFD](https://github.com/warwick-hpsc/MG-CFD-app-plain). Provides MPI, full OpenMP, SIMD, CUDA, OpenACC, OpenMP 4, and some pairings of them.

Dependencies
==========================================

* [OP2](https://github.com/OP-DSL/OP2-Common)
* HDF5 parallel
* ParMETIS
* PT-Scotch
* CUDA library if using `mgcfd_cuda` or `mgcfd_mpi_cuda`

Compiling and executing
==========================================

### Building dependencies:

#### OP2
Several distinct libraries can be compiled, depending on the mix of parallelism and performance portability desired. Similary for MG-CFD, and so it is likely that you only need to compile a subset of the OP2 libraries. 

All MG-CFD variants need two particular OP2 libraries, created by executing these two OP2 make rules: `core`, `hdf5`. Variant-specific dependencies are listed in next section.

#### HDF5
All variants of MG-CFD require HDF5. Buy only those MPI variants require HDF5 compiled with MPI support (with --enable-parallel)

#### MPI-specific
MPI variants of MG-CFD require these further libraries

##### ParMETIS
Standard compile.

##### PT-Scotch
Follow their build instructions. After linking `Makefile.inc`, edit it and remove the flag `-DSCOTCH_PTHREAD` from `CFLAGS`. Then standard compile.

### Compiling MG-CFD

Different binaries can be generated, depending on the mix of parallelism and performance portability desired:

Intent | MG-CFD make rule | OP2 dependency make rule
------ | --------- | -----------------------------
Sequential | seq | seq
OpenMP | openmp | openmp
MPI | mpi | mpi_seq
MPI + OpenMP | mpi_openmp | mpi_seq
MPI + SIMD | mpi_vec | mpi_seq
CUDA | cuda | cuda
MPI + CUDA | mpi_cuda | mpi_cuda

In future, OpenACC and OpenMP 4.5 ports will be available

### Quick run:

Want to execute immediately? Navigate to a folder containing input HDF5 files and execute:

```Shell
     $ ./path/to/mgcfd_* -i input.dat
```

MG-CFD has more command-line arguments to ease file/directory interaction, and control execution parameters. View the help page for more information:

```Shell
     $ ./path/to/mgcfd_* --help
```

### Performance counters:

Built into MG-CFD is functionality to collect performance counter data, at fine granularity of individual loops. Currently CPU only. Requires PAPI library to be installed and configured. 
- disabled as default - to enable, enable either 'PAPI' flag in Makefile, then compile. 
- this in turn enables a command-line parameter: -p \<filepath\> . This file should contain the list of events to measure.
- counts will be written to PAPI.csv

### Generating batch submission scripts:

1) Prepare a json file detailing run configuration. See ./run-inputs/annotated.json for documentation on each option. 

2) Generate run batch scripts from the json file:

```Shell
     $ python ./run-scripts/gen_job.py --json path/to/config.json
```
     
3) The specified `jobs directory` will contain a subfolder for each run configuration, and a single `submit_all.sh` file. If a scheduler was specified in the .json file, then `submit_all.sh` will compile locally then submit each job to the scheduler for execution. If local execution was requested in the json file, then `submit_all.sh` will compile and execute locally. 

4) Each run will output .csv files containing performance data. These can be collated together using `aggregate-output-data.py`:

```Shell
     $ python ./run-scripts/aggregate-output-data.py \ 
              --output-dirpath path/to/desired-folder/for/collated-csv-files \
              --data-dirpath path/to/job-output
```

Validating result
==========================================

MG-CFD can verify the final flow state against a precomputed solution file, useful for assuring correctness of code changes. To perform this use the `-v` parameter, and set the number of multigrid cycles `-g` to match the solution file (inspect its filename).

```Shell
     $ mgcfd_* ... -v -g 10
```

You can also generate your own solution file:

```Shell
     $ mgcfd_* ... --output-variables --output-file-prefix "solution."
```

This will generate a solution file for each multigrid level, e.g. `solution.variables.L0.cycles=10.h5`

Datasets
==========================================

A [release](https://github.com/warwick-hpsc/MG-CFD-app-OP2/releases) is provided that includes the [Onera M6 wing](https://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html). It consists of 300K nodes (930K edges), and three additional multigrid meshes with respective node counts of 165K, 111K, and 81K. 

Additional larger meshes are available at [our research groups's homepage](https://warwick.ac.uk/fac/sci/dcs/research/systems/hpsc/software):
* Rotor 37 1M cells (multigrid)
* Rotor 37 8M cells (multigrid)
* Rotor 37 25M cells (multigrid)
* Rotor 37 150M cells (single level)

Major updates since release
==========================================
12/Jun/2019: added MPI + SIMD variant

Authorship
==========================================

Andrew Owenson: a.owenson@warwick.ac.uk

For more information on the design of MG-CFD, please refer to our publication: https://onlinelibrary.wiley.com/doi/10.1002/cpe.5443

If you wish to cite this work then please use the following:

* Owenson A.M.B., Wright S.A., Bunt R.A., Ho Y.K., Street M.J., and Jarvis S.A. (2019), An Unstructured CFD Mini-Application for the Performance Prediction of a Production CFD Code, Concurrency Computat: Pract Exper., 2019

