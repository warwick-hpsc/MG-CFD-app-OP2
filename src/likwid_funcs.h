//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//


#ifndef LIKWID_FUNCS_H
#define LIKWID_FUNCS_H

#include <likwid.h>
#include <vector>

#define long_long long long

// OP2:
#include "op_seq.h"

#include "config.h"
#include "const.h"
#include "global.h"

extern int n_events;
extern char** event_names;
extern int likwid_gid;
extern CpuTopology_t likwid_topo;

extern long_long** flux_kernel_event_counts;
extern long_long** ustream_kernel_event_counts;

inline int likwid_get_num_threads() {
    #ifdef _OMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

inline void my_likwid_start()
{
	if (n_events == 0) {
		return;
	}
	int err = perfmon_startCounters();
	if (err < 0) {
		op_printf("Failed to start counters for group %d for thread %d\n", likwid_gid, (-1*err)-1);
		exit(EXIT_FAILURE);
	}
}

inline void my_likwid_stop(long_long** restrict event_counts)
{
	if (n_events == 0) {
		return;
	}
	int err = perfmon_stopCounters();
	if (err < 0) {
		op_printf("Failed to stop counters for group %d for thread %d\n",likwid_gid, (-1*err)-1);
		exit(EXIT_FAILURE);
	}
	for (int e=0; e<n_events; e++) {
		for (int t=0; t<likwid_topo->numHWThreads; t++) {
			double result_f = perfmon_getResult(likwid_gid, e, t);
			long_long result = (long_long)result_f;
			event_counts[t][current_level*n_events + e] += result;
		}
	}
}

inline void init_likwid()
{
    int err, j;
    err = topology_init();
    if (err < 0) {
        op_printf("Failed to initialize LIKWID's topology module\n");
        exit(EXIT_FAILURE);
    }
    CpuInfo_t info = get_cpuInfo();
    likwid_topo = get_cpuTopology();
    affinity_init();

    int* cpus = (int*)malloc(likwid_topo->numHWThreads * sizeof(int));
    if (!cpus) {
        exit(EXIT_FAILURE);
    }
    for (j=0; j<likwid_topo->numHWThreads;j++) {
        cpus[j] = likwid_topo->threadPool[j].apicId;
    }

    err = perfmon_init(likwid_topo->numHWThreads, cpus);
    if (err < 0) {
        op_printf("Failed to initialize LIKWID's performance monitoring module\n");
        topology_finalize();
        exit(EXIT_FAILURE);
    }

    // numa_init();
    // perfmon_init_maps();
    // perfmon_check_counter_map(0);

    n_events = 0;
    likwid_gid = 0;
}

inline void load_likwid_events()
{
    std::vector<std::string> events;

    std::string line;
    std::ifstream file_if(conf.papi_config_file);
    if(!file_if.is_open()) {
        op_printf("WARNING: Failed to open Likwid config: '%s'\n", conf.papi_config_file);
        exit(EXIT_FAILURE);
    } else {
        while(std::getline(file_if, line)) {
            if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
                continue;
            }
            events.push_back(line);
        }
    }

    n_events = events.size();
    if (n_events > 0) {
    	int j, err, n;

    	// Store event names
    	event_names = (char**)malloc(n_events*sizeof(char*));
    	for (j=0; j<n_events; j++) {
    		n = strlen(events[j].c_str()) + 1;
    		event_names[j] = (char*)malloc(n*sizeof(char));
    		strcpy(event_names[j], events[j].c_str());
    	}

    	// Construct comma-delimited event string
    	int estr_length = strlen(events[0].c_str());
		for (j=1; j<n_events; j++) {
			estr_length += 1 + strlen(events[j].c_str());
		}
		char* estr = (char*)malloc(estr_length*sizeof(char));
		strcpy(estr, events[0].c_str());
		int i = 0;
		for (j=1; j<n_events; j++) {
			i += strlen(events[j-1].c_str());
			estr[i] = ',';
			i++;
			strcpy(estr+i, events[j].c_str());
		}

		// Setup Likwid
		likwid_gid = perfmon_addEventSet(estr);
		if (likwid_gid < 0) {
			op_printf("Failed to add event string %s to LIKWID's performance monitoring module\n", estr);
			perfmon_finalize();
			topology_finalize();
			exit(EXIT_FAILURE);
		}
		err = perfmon_setupCounters(likwid_gid);
		if (err < 0) {
			op_printf("Failed to setup group %d in LIKWID's performance monitoring module\n", likwid_gid);
			perfmon_finalize();
			topology_finalize();
			exit(EXIT_FAILURE);
		}

	    const int nt = likwid_topo->numHWThreads;
	    flux_kernel_event_counts = (long_long**)malloc(sizeof(long_long*)*nt);
    	ustream_kernel_event_counts = (long_long**)malloc(sizeof(long_long*)*nt);
    	for (int tid=0; tid<nt; tid++) {
	        flux_kernel_event_counts[tid] = (long_long*)malloc(sizeof(long_long)*levels*n_events);
	        for (int i=0; i<levels*n_events; i++) {
	            flux_kernel_event_counts[tid][i] = 0;
	        }
	        ustream_kernel_event_counts[tid] = (long_long*)malloc(sizeof(long_long)*levels*n_events);
	        for (int i=0; i<levels*n_events; i++) {
	            ustream_kernel_event_counts[tid][i] = 0;
	        }
    	}

		op_printf("LIKWID initialised, monitoring '%s'\n", estr);
    }
}

inline void dump_likwid_counters_to_file(
	int rank, 
    long_long** flux_kernel_event_counts, 
    long_long** ustream_kernel_event_counts, 
    char* output_file_prefix)
{
	std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += std::string("P=") + number_to_string(rank);
    filepath += ".Likwid.csv";

    bool write_header = false;
    std::ifstream f(filepath.c_str());
    if (!f || f.peek() == std::ifstream::traits_type::eof()) {
        write_header = true;
    }
    f.close();

    std::ostringstream header;

    const int nt = likwid_get_num_threads();
    if (write_header) {
        header << "Rank";
        header << ",Thread";
        header << ",Partitioner";
        header << ",Likwid event";
        header << ",kernel";
        header << ",level";
        header << ",count";
    }

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::app);
    if (write_header) outfile << header.str() << std::endl;

    for (int tid=0; tid<nt; tid++) {
        for (int eid=0; eid<n_events; eid++) {
            for (int l=0; l<levels; l++) {
                std::ostringstream event_data_line;
                event_data_line << rank;
                event_data_line << "," << tid;
                event_data_line << "," << conf.partitioner_string;
                event_data_line << "," << event_names[eid];

                event_data_line << "," << "compute_flux_edge_kernel";
                event_data_line << "," << l;

                const int idx = l*n_events + eid;
                event_data_line << ',' << flux_kernel_event_counts[tid][idx];

                outfile << event_data_line.str() << std::endl;
            }

            if (conf.measure_mem_bound) {
                for (int l=0; l<levels; l++) {
                    std::ostringstream event_data_line;
                    event_data_line << rank;
                    event_data_line << "," << tid;
                    event_data_line << "," << conf.partitioner_string;
                    event_data_line << "," << event_names[eid];

                    event_data_line << "," << "unstructured_stream_kernel";
                    event_data_line << "," << l;

                    const int idx = l*n_events + eid;
                    event_data_line << ',' << ustream_kernel_event_counts[tid][idx];

                    outfile << event_data_line.str() << std::endl;
                }
            }
        }
    }
    outfile.close();
}

#endif