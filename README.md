MG-CFD OP2
==========================================

OP2 port of [MG-CFD](https://github.com/warwick-hpsc/MG-CFD-app-plain). Includes all accelerations provided by OP2.

Dependencies
==========================================

* [OP2](https://github.com/OP-DSL/OP2-Common) including HDF5 support.

Compiling and executing
==========================================

### Building dependencies:

#### OP2

Follow the OP2 documentation to install the OP2 library. Building the HDF5 variants of OP2 is required.

Several backends can be compiled, depending on the desired combination of parallelism and performance portability.

### Compiling MG-CFD

Different binaries can be generated depending on the desired combination of parallelism and performance portability.

Simply run `make` within the repository folder to compile all available MG-CFD variants.

Note: You must set the environment variable `OP2_INSTALL_PATH` to point to the OP2 installation directory.

If you have just compiled OP2, you can set it as follows:

```
export OP2_INSTALL_PATH=<path-to-op2-common>/op2
```

### Quick run:

To run immediately, navigate to a directory containing input HDF5 files and execute:

```Shell
     $ ./<mgcfd-repo-path>/euler3d_* -i input.dat
```

MG-CFD provides additional command-line arguments for file/directory handling and execution control. See the help page for more details:

```Shell
     $ ./<mgcfd-repo-path>/euler3d_* --help
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
     $ <mgcfd-repo-path>/euler3d_* ... -v -g 10
```

You can also generate your own solution file:

```Shell
     $ <mgcfd-repo-path>/euler3d_* ... --output-variables --output-file-prefix "solution."
```

This will generate a solution file for each multigrid level, e.g. `solution.variables.L0.cycles=10.h5`

Datasets
==========================================

A [release](https://github.com/warwick-hpsc/MG-CFD-app-OP2/releases) is provided that includes the [Onera M6 wing](https://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html). It consists of 300K nodes (930K edges), and three additional multigrid meshes with respective node counts of 165K, 111K, and 81K. 

Additional larger meshes are available at [our research group's homepage](https://warwick.ac.uk/fac/sci/dcs/research/systems/hpsc/software):
* Rotor 37 1M cells (multigrid)
* Rotor 37 8M cells (multigrid)
* Rotor 37 25M cells (multigrid)
* Rotor 37 150M cells (single level)

Major updates since release
==========================================

12/Jun/2019: added MPI + SIMD variant

18/March/2026: Support OP2 new code-generator

Authorship
==========================================

Andrew Owenson: a.owenson@warwick.ac.uk

For more information on the design of MG-CFD, please refer to our publication: https://onlinelibrary.wiley.com/doi/10.1002/cpe.5443

If you wish to cite this work then please use the following:

* Owenson A.M.B., Wright S.A., Bunt R.A., Ho Y.K., Street M.J., and Jarvis S.A. (2019), An Unstructured CFD Mini-Application for the Performance Prediction of a Production CFD Code, Concurrency Computat: Pract Exper., 2019

