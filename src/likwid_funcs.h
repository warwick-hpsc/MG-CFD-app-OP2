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

#define long_long long long

// OP2:
#include "op_seq.h"

#include "config.h"

extern int n_events;
extern int likwid_gid;

inline int likwid_get_num_threads() {
    #ifdef _OMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

inline void init_likwid()
{
        int err, j;
        err = topology_init();
        if (err < 0) {
            printf("Failed to initialize LIKWID's topology module\n");
            return;
        }
        CpuInfo_t info = get_cpuInfo();
        CpuTopology_t topo = get_cpuTopology();
        affinity_init();

        int* cpus = (int*)malloc(topo->numHWThreads * sizeof(int));
        if (!cpus) {
            return;
        }
        for (j=0; j<topo->numHWThreads;j++) {
            cpus[j] = topo->threadPool[j].apicId;
        }

        err = perfmon_init(topo->numHWThreads, cpus);
        if (err < 0) {
            printf("Failed to initialize LIKWID's performance monitoring module\n");
            topology_finalize();
            return;
        }

        // numa_init();
        // perfmon_init_maps();
        // perfmon_check_counter_map(0);

        n_events = 0;
        likwid_gid = 0;
}

inline void load_likwid_events()
{
    // num_events = (int*)malloc(sizeof(int)*likwid_get_num_threads());

	std::vector<std::string> events;

    std::string line;
    std::ifstream file_if(conf.papi_config_file);
    if(!file_if.is_open()) {
        op_printf("WARNING: Failed to open Likwid config: '%s'\n", conf.papi_config_file);
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
    	int estr_length = strlen(events[0].c_str());
    	int j, err;
		for (j=1; j<n_events; j++) {
			estr_length += 1 + strlen(events[j].c_str());
		}

		char* estr = NULL;
		estr = (char*)malloc(estr_length*sizeof(char));
		strcpy(estr, events[0].c_str());
		int i = 0;
		for (j=1; j<n_events; j++) {
			i += strlen(events[j-1].c_str());
			estr[i] = ',';
			i++;
			strcpy(estr+i, events[j].c_str());
		}

		likwid_gid = perfmon_addEventSet(estr);
		if (likwid_gid < 0) {
			printf("Failed to add event string %s to LIKWID's performance monitoring module\n", estr);
			perfmon_finalize();
			topology_finalize();
			return;
		}
		err = perfmon_setupCounters(likwid_gid);
		if (err < 0) {
			printf("Failed to setup group %d in LIKWID's performance monitoring module\n", likwid_gid);
			perfmon_finalize();
			topology_finalize();
			return;
		}

		printf("LIKWID initialised, monitoring '%s'\n", estr);
    }
}

#endif