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


#ifndef PAPI_H
#define PAPI_H

#ifdef PAPI

#include <string>
#include <fstream>
#include <sstream>

#include <papi.h>

// OP2:
#include "op_seq.h"

#include "config.h"

inline void my_papi_start(int event_set)
{
    if (event_set != PAPI_NULL) {
        if (PAPI_start(event_set) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to start PAPI\n");
            exit(EXIT_FAILURE);
        }
    }
}

inline void my_papi_stop(
    long_long* restrict event_counts, 
    long_long* restrict temp_count_store, 
    int event_set, 
    int num_events)
{
	if (event_set != PAPI_NULL) {
	    if (PAPI_stop(event_set, temp_count_store) != PAPI_OK) {
	    	fprintf(stderr, "ERROR: Failed to stop PAPI\n");
	        exit(EXIT_FAILURE);
	    }
	    for (int e=0; e<num_events; e++) {
	        event_counts[e] += temp_count_store[e];
	    }
	}
}

inline unsigned long omp_get_thread_num_ul() {
    #ifdef _OMP
        return (unsigned long)omp_get_thread_num();
    #else
        return 0;
    #endif
}

inline void init_papi(int* num_events)
{
    int ret;

    if (PAPI_num_counters() < 2) {
       fprintf(stderr, "ERROR: No hardware counters here, or PAPI not supported.\n");
       exit(-1);
    }

    std::string line;
    std::ifstream file_if(conf.papi_config_file);
    *num_events = 0;
    while(std::getline(file_if, line))
    {
        if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
            continue;
        }
        (*num_events)++;
    }

    if ((ret=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
        fprintf(stderr, "ERROR: PAPI_library_init() failed: '%s'.\n", PAPI_strerror(ret));
        exit(EXIT_FAILURE);
    }

    if ((ret=PAPI_thread_init(omp_get_thread_num_ul)) != PAPI_OK) {
        fprintf(stderr, "ERROR: PAPI_thread_init() failed: '%s'.\n", PAPI_strerror(ret));
        exit(EXIT_FAILURE);
    }
}

inline void load_papi_events(int num_events, int* event_set, int** events)
{
    int ret;

    *event_set = PAPI_NULL;
    ret = PAPI_create_eventset(event_set);
    if (ret != PAPI_OK || (*event_set)==PAPI_NULL) {
        fprintf(stderr, "ERROR: PAPI_create_eventset() failed: '%s'.\n", PAPI_strerror(ret));
        exit(EXIT_FAILURE);
    }

    char event_name[512];
    std::string line;
    std::ifstream file_if(conf.papi_config_file);
    if(!file_if.is_open()) {
        // op_printf("ERROR: Failed to open PAPI config: '%s'\n", conf.papi_config_file);
        // exit(EXIT_FAILURE);
        op_printf("WARNING: Failed to open PAPI config: '%s'\n", conf.papi_config_file);
        *event_set = PAPI_NULL;
        return;
    }
    while(std::getline(file_if, line))
    {
        if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
            continue;
        }

        strcpy(event_name, line.c_str());
        // printf("Processing event: %s\n", event_name);

        // ret = PAPI_add_named_event((*event_set), event_name);
        ret = PAPI_OK - 1;
        if (ret != PAPI_OK) {
            // printf("PAPI_add_named_event() failed, attempting add by code conversion\n"); fflush(stdout);

            // fprintf(stderr, "Thread %d on CPU %d: failed to add event %s to event set: '%s'.\n", omp_get_thread_num(), sched_getcpu(), event_name, PAPI_strerror(ret));
            // if ((*event_set)==PAPI_NULL) {
            //     fprintf(stderr, "... and event_set=PAPI_NULL\n");
            // }

            // It appears that PAPI_add_named_event() only works with native events, not PAPI presets.
            int code = -1;
            ret = PAPI_event_name_to_code(event_name, &code);
            if (ret != PAPI_OK) {
                op_printf("Could not convert string '%s' to PAPI event, error = %s\n", event_name, PAPI_strerror(ret));
            } else {
                if (PAPI_query_event(code) != PAPI_OK) {
                    op_printf("PAPI event %s not present\n", event_name);
                } else {
                    ret = PAPI_add_event((*event_set), code);
                    if (ret != PAPI_OK) {
                        fprintf(stderr, "ERROR: Failed to add event %d to event set: '%s'.\n", code, PAPI_strerror(ret));
                        if ((*event_set)==PAPI_NULL) {
                            fprintf(stderr, "ERROR: ... and event_set=PAPI_NULL\n");
                        }
                        exit(EXIT_FAILURE);
                    }
                    // else {
                    //     printf("Monitoring PAPI event '%s'\n", event_name);
                    // }
                }
            }
        }
    }
    if (file_if.bad()) {
        op_printf("ERROR: Failed to read papi conf file: %s\n", conf.papi_config_file);
        exit(EXIT_FAILURE);
    }
    // printf("Finished parsing PAPI file\n");

    *events = (int*)malloc(sizeof(int)*num_events);
    if (PAPI_list_events(*event_set, *events, &num_events) != PAPI_OK) {
        fprintf(stderr, "ERROR: PAPI_list_events() failed\n");
        exit(EXIT_FAILURE);
    }

    // long_long* temp_count_stores = alloc<long_long>(num_events);
    long_long* temp_count_stores = (long_long*)malloc(sizeof(long_long)*num_events);
    for (int e=0; e<num_events; e++) {
        temp_count_stores[e] = 0;
    }
    if ((*event_set) != PAPI_NULL) {
        if (PAPI_start((*event_set)) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to start PAPI\n");
            exit(EXIT_FAILURE);
        }
        if (PAPI_stop((*event_set), temp_count_stores) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to stop PAPI\n");
            exit(EXIT_FAILURE);
        }
    }
    free(temp_count_stores);
}

inline void dump_papi_counters_to_file(
	int rank, 
    int num_levels, 
    int num_events,
    int* events,  
    long_long* flux_kernel_event_counts, 
    char* output_file_prefix)
{
    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += std::string("P=") + number_to_string(rank);
    filepath += ".PAPI.csv";

    bool write_header = false;
    std::ifstream f(filepath.c_str());
    if (!f || f.peek() == std::ifstream::traits_type::eof()) {
        write_header = true;
    }
    f.close();

    std::ostringstream header;

    if (write_header) {
        header << "Rank";
        header << ",Partitioner";
        header << ",PAPI counter";
        header << ",kernel";
        header << ",level";
        header << ",count";
    }

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::app);
    if (write_header) outfile << header.str() << std::endl;

    for (int eid=0; eid<num_events; eid++)
    {
        char eventName[PAPI_MAX_STR_LEN] = "";
        if (PAPI_event_code_to_name(events[eid], eventName) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to convert code %d to name\n", events[eid]);
            exit(EXIT_FAILURE);
        }

        for (int l=0; l<num_levels; l++) {
            std::ostringstream event_data_line;
            event_data_line << rank;
            event_data_line << "," << conf.partitioner_string;
            event_data_line << "," << eventName;

            event_data_line << "," << "compute_flux_edge_kernel";
            event_data_line << "," << l;

            const int idx = l*num_events + eid;
            event_data_line << ',' << flux_kernel_event_counts[idx];

            outfile << event_data_line.str() << std::endl;
        }
    }
    outfile.close();
}

#endif

#endif
