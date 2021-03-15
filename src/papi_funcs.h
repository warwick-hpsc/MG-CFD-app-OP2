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

// OP2 op_printf
#include "op_seq.h"

#include "config.h"
#include "const.h"
#include "global.h"

extern int* num_events;
extern int* event_set;
extern int** events;
extern long_long** flux_kernel_event_counts;
extern long_long** ustream_kernel_event_counts;
extern long_long** temp_count_store;

typedef struct papi_row_t
{
    int thread_id;
    int cpu_id;
    char event_name[TABLE_STRING_LENGTH];
    char kernel_name[TABLE_STRING_LENGTH];
    int level;
    long long int count;
} papi_row;

#ifdef _MPI
#include <mpi.h>
inline void dump_papi_counters_to_file_mpi(
    papi_row* my_table, int n_rows, 
    char* output_file_prefix);
#endif

inline unsigned long omp_get_thread_num_ul() {
    #ifdef _OMP
        return (unsigned long)omp_get_thread_num();
    #else
        return 0;
    #endif
}

inline int papi_get_num_threads() {
    #ifdef _OMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

inline int papi_get_thread_id() {
    #ifdef _OMP
        return omp_get_thread_num();
    #else
        return 0;
    #endif
}

inline void my_papi_start()
{
    const int tid = papi_get_thread_id();
    if (event_set[tid] != PAPI_NULL) {
        for (int e=0; e<num_events[tid]; e++) {
            temp_count_store[tid][e] = 0;
        }
        if (PAPI_start(event_set[tid]) != PAPI_OK) {
            fprintf(stderr, "ERROR: Failed to start PAPI\n");
            exit(EXIT_FAILURE);
        }
    }
}

inline void my_papi_stop(long_long** restrict event_counts)
{
    const int tid = papi_get_thread_id();
	if (event_set[tid] != PAPI_NULL) {
	    if (PAPI_stop(event_set[tid], temp_count_store[tid]) != PAPI_OK) {
	    	fprintf(stderr, "ERROR: Failed to stop PAPI\n");
	        exit(EXIT_FAILURE);
	    }
        const int n = num_events[tid];
	    for (int e=0; e<n; e++) {
	        event_counts[tid][current_level*n + e] += temp_count_store[tid][e];
	    }
	}
}

inline void init_papi()
{
    int ret;

    if ((ret=PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
        fprintf(stderr, "ERROR: PAPI_library_init() failed: '%s'.\n", PAPI_strerror(ret));
        exit(EXIT_FAILURE);
    }

    // if (PAPI_num_counters() < 2) {
    int num_ctrs = 0;
    int num_comps = PAPI_num_components();
    for (int c=0; c<num_comps; c++) {
        num_ctrs += PAPI_num_cmp_hwctrs(c);
    }
    if (num_ctrs < 2) {
       fprintf(stderr, "ERROR: No hardware counters here, or PAPI not supported (num_ctrs=%d)\n", num_ctrs);
       exit(-1);
    }

    #ifdef _OMP
        #pragma omp parallel
        {
            #pragma omp critical
            {
                if ((ret=PAPI_thread_init(omp_get_thread_num_ul)) != PAPI_OK) {
                    fprintf(stderr, "ERROR: PAPI_thread_init() failed: '%s'.\n", PAPI_strerror(ret));
                    exit(EXIT_FAILURE);
                }
            }
        }
    #else
        if ((ret=PAPI_thread_init(omp_get_thread_num_ul)) != PAPI_OK) {
            fprintf(stderr, "ERROR: PAPI_thread_init() failed: '%s'.\n", PAPI_strerror(ret));
            exit(EXIT_FAILURE);
        }
    #endif
}

inline void load_papi_events()
{
    int ret;

    event_set = (int*)malloc(sizeof(int)*papi_get_num_threads());
    num_events = (int*)malloc(sizeof(int)*papi_get_num_threads());

    bool first_thread_is_initialising = true;
    #ifdef _OMP
        #pragma omp parallel private(ret)
        {
            const int tid = omp_get_thread_num();
            #pragma omp critical
            {
    #else
        const int tid = 0;
    #endif

        event_set[tid] = PAPI_NULL;
        ret = PAPI_create_eventset(&event_set[tid]);
        if (ret != PAPI_OK || event_set[tid]==PAPI_NULL) {
            fprintf(stderr, "ERROR: PAPI_create_eventset() failed: '%s'.\n", PAPI_strerror(ret));
            exit(EXIT_FAILURE);
        }

        char event_name[512];
        std::string line;
        std::ifstream file_if(conf.papi_config_file);
        if(!file_if.is_open()) {
            op_printf("WARNING: Failed to open PAPI config: '%s'\n", conf.papi_config_file);
            event_set[tid] = PAPI_NULL;
        } else { 
            while(std::getline(file_if, line)) {
                if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
                    continue;
                }

                if (line.find(std::string("unc_imc")) != std::string::npos) {
                    // This is a uncore imc event. If many threads monitor it 
                    // then PAPI throws up errors, so allow only one thread to 
                    // monitor it.
                    if (!first_thread_is_initialising) {
                        continue;
                    }
                }
                if (line.find(std::string("unc_edc")) != std::string::npos) {
                    // This is a uncore mcdram event. If many threads monitor it 
                    // then PAPI throws up errors, so allow only one thread to 
                    // monitor it.
                    if (!first_thread_is_initialising) {
                      continue;
                    }
                }

                strcpy(event_name, line.c_str());
                // printf("Processing event: %s\n", event_name);

                // ret = PAPI_add_named_event((event_set[tid]), event_name);
                ret = PAPI_OK - 1;
                if (ret != PAPI_OK) {
                    // printf("PAPI_add_named_event() failed, attempting add by code conversion\n"); fflush(stdout);

                    // fprintf(stderr, "Thread %d on CPU %d: failed to add event %s to event set: '%s'.\n", omp_get_thread_num(), sched_getcpu(), event_name, PAPI_strerror(ret));
                    // if (event_set[tid]==PAPI_NULL) {
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
                            ret = PAPI_add_event(event_set[tid], code);
                            if (ret != PAPI_OK) {
                                fprintf(stderr, "ERROR: Failed to add event %d to event set: '%s'.\n", code, PAPI_strerror(ret));
                                if (event_set[tid]==PAPI_NULL) {
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
        }
        if (file_if.bad()) {
            op_printf("ERROR: Failed to read papi conf file: %s\n", conf.papi_config_file);
            exit(EXIT_FAILURE);
        }
        // printf("Finished parsing PAPI file\n");

        int n = PAPI_num_events(event_set[tid]);
        num_events[tid] = n;
        if (n == 0) {
            event_set[tid] = PAPI_NULL;
        }
        
        first_thread_is_initialising = false;

    #ifdef _OMP
            } // End of critical region
        } // End of parallel loop
    #endif

    const int nt = papi_get_num_threads();
    events = (int**)malloc(sizeof(int*)*nt);
    flux_kernel_event_counts = (long_long**)malloc(sizeof(long_long*)*nt);
    ustream_kernel_event_counts = (long_long**)malloc(sizeof(long_long*)*nt);
    temp_count_store = (long_long**)malloc(sizeof(long_long*)*nt);
    for (int tid=0; tid<papi_get_num_threads(); tid++) {
        int n = num_events[tid];

        events[tid] = (int*)malloc(sizeof(int)*n);
        int* temp_thread_events = (int*)malloc(sizeof(int)*n);
        if (PAPI_list_events(event_set[tid], temp_thread_events, &n) != PAPI_OK) {
            fprintf(stderr, "ERROR: PAPI_list_events() failed\n");
            DEBUGGABLE_ABORT
        }
        for (int i=0; i<n; i++) {
            events[tid][i] = temp_thread_events[i];
        }
        free(temp_thread_events);

        flux_kernel_event_counts[tid] = (long_long*)malloc(sizeof(long_long)*levels*n);
        for (int i=0; i<levels*n; i++) {
            flux_kernel_event_counts[tid][i] = 0;
        }

        ustream_kernel_event_counts[tid] = (long_long*)malloc(sizeof(long_long)*levels*n);
        for (int i=0; i<levels*n; i++) {
            ustream_kernel_event_counts[tid][i] = 0;
        }

        temp_count_store[tid] = (long_long*)malloc(sizeof(long_long)*n);
        for (int e=0; e<n; e++) {
            temp_count_store[tid][e] = 0;
        }
    }

    #ifdef _OMP
        #pragma omp parallel
        {
            const int tid = papi_get_thread_id();
            int n = num_events[tid];
    #endif
        if (event_set[tid] != PAPI_NULL) {
            if (PAPI_start(event_set[tid]) != PAPI_OK) {
                fprintf(stderr, "ERROR: Thread %d failed to start PAPI\n", tid);
                exit(EXIT_FAILURE);
            }
            if (PAPI_stop(event_set[tid], temp_count_store[tid]) != PAPI_OK) {
                fprintf(stderr, "ERROR: Thread %d failed to stop PAPI\n", tid);
                exit(EXIT_FAILURE);
            }
            for (int e=0; e<n; e++) {
                temp_count_store[tid][e] = 0;
            }
        }
    #ifdef _OMP
        } // End of parallel loop
    #endif
}

inline void dump_papi_counters_to_file(
    long_long** flux_kernel_event_counts, 
    long_long** ustream_kernel_event_counts, 
    char* output_file_prefix)
{
    const int nt = papi_get_num_threads();
    int n_events = 0;
    for (int tid=0; tid<nt; tid++) {
        n_events += num_events[tid];
    }
    int n_rows = nt * n_events * levels;
    if (conf.measure_mem_bound) {
        n_rows *= 2;
    }
    papi_row* my_table = (papi_row*)malloc(n_rows*sizeof(papi_row));
    int* cpu_ids = (int*)malloc(nt*sizeof(int));
    #ifdef _OMP
        #pragma omp parallel
        {
            cpu_ids[omp_get_thread_num()] = sched_getcpu();
        }
    #else
        cpu_ids[0] = sched_getcpu();
    #endif

    int row_idx = 0;
    for (int tid=0; tid<nt; tid++) {
        for (int eid=0; eid<num_events[tid]; eid++) {
            char eventName[PAPI_MAX_STR_LEN] = "";
            if (PAPI_event_code_to_name(events[tid][eid], eventName) != PAPI_OK) {
                fprintf(stderr, "ERROR: Failed to convert code %d to name\n", events[tid][eid]);
                exit(EXIT_FAILURE);
            }

            for (int l=0; l<levels; l++) {
                my_table[row_idx].thread_id = tid;
                my_table[row_idx].cpu_id = cpu_ids[tid];
                strncpy(my_table[row_idx].event_name, eventName, TABLE_STRING_LENGTH);
                strncpy(my_table[row_idx].kernel_name, "compute_flux_edge_kernel", TABLE_STRING_LENGTH);
                my_table[row_idx].level = l;

                const int idx = l*num_events[tid] + eid;
                my_table[row_idx].count = flux_kernel_event_counts[tid][idx];
                row_idx++;
            }

            if (conf.measure_mem_bound) {
                for (int l=0; l<levels; l++) {
                    my_table[row_idx].thread_id = tid;
                    my_table[row_idx].cpu_id = cpu_ids[tid];
                    strncpy(my_table[row_idx].event_name, eventName, TABLE_STRING_LENGTH);
                    strncpy(my_table[row_idx].kernel_name, "unstructured_stream_kernel", TABLE_STRING_LENGTH);
                    my_table[row_idx].level = l;

                    const int idx = l*num_events[tid] + eid;
                    my_table[row_idx].count = ustream_kernel_event_counts[tid][idx];
                    row_idx++;
                }
            }
        }
    }

    #ifdef _MPI
        dump_papi_counters_to_file_mpi(my_table, n_rows, output_file_prefix);
        return;
    #endif

    // Write out data:
    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "PAPI.csv";

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::out);
    std::ostringstream header;
    if (nt > 1) {
        header << ",Thread";
        header << ",CpuId";
    }
    header << ",Partitioner";
    header << ",PAPI counter";
    header << ",kernel";
    header << ",level";
    header << ",count";
    outfile << header.str() << std::endl;

    for (int row_idx=0; row_idx<n_rows; row_idx++) {
        std::ostringstream event_data_line;
        if (nt > 1) {
            event_data_line << "," << my_table[row_idx].thread_id;
            event_data_line << "," << my_table[row_idx].cpu_id;
        }
        event_data_line << "," << conf.partitioner_string;
        event_data_line << "," << my_table[row_idx].event_name;
        event_data_line << "," << my_table[row_idx].kernel_name;
        event_data_line << "," << my_table[row_idx].level;
        event_data_line << "," << my_table[row_idx].count;
        outfile << event_data_line.str() << std::endl;
    }

    outfile.close();
}

// Functionality below should be moved into OP2
#ifdef _MPI
#include <mpi.h>
inline void dump_papi_counters_to_file_mpi(
    papi_row* my_table, int n_rows, 
    char* output_file_prefix)
{
    // Non-root ranks send data to root, which writes to file.

    int err;

    // Create MPI_Datatype MPI_PapiRowType
    int types[6] = { MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_LONG_LONG_INT};
    int blocklengths[6] = { 1, 1, TABLE_STRING_LENGTH, TABLE_STRING_LENGTH, 1, 1 };
    int type_sizes[6] = { sizeof(int), sizeof(int), sizeof(char), sizeof(char), sizeof(int), sizeof(long long int) };
    MPI_Aint displacements[6];
    displacements[0] = offsetof(papi_row, thread_id);
    displacements[1] = offsetof(papi_row, cpu_id);
    displacements[2] = offsetof(papi_row, event_name);
    displacements[3] = offsetof(papi_row, kernel_name);
    displacements[4] = offsetof(papi_row, level);
    displacements[5] = offsetof(papi_row, count);
    MPI_Datatype MPI_PapiRowType;
    err = MPI_Type_create_struct(
        6,
        blocklengths,
        displacements,
        types,
        &MPI_PapiRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to create custom MPI 'PAPI table' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }
    err = MPI_Type_commit(&MPI_PapiRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to commit custom MPI 'PAPI row' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    // Send/receive tables
    int nranks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    papi_row** papi_tables = NULL;
    int* n_rows_each = (int*)malloc(nranks * sizeof(int));
    err = MPI_Gather(&n_rows, 1, MPI_INT, n_rows_each, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: MPI_Irecv() failed\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        papi_tables = (papi_row**)malloc(nranks*sizeof(papi_row*));
        for (int r=0; r<nranks; r++) {
            papi_tables[r] = (papi_row*)malloc(n_rows_each[r]*sizeof(papi_row));
        }

        MPI_Request* reqs = (MPI_Request*)malloc((nranks-1)*sizeof(MPI_Request));

        for (int r=1; r<nranks; r++) {
            err = MPI_Irecv(papi_tables[r], n_rows_each[r], MPI_PapiRowType, r, 0, MPI_COMM_WORLD, &reqs[r-1]);
            if (err != MPI_SUCCESS) {
                op_printf("ERROR: MPI_Irecv() failed\n");
                op_exit();
                exit(EXIT_FAILURE);
            }
        }
        int num_reqs_pending = (nranks-1);
        while (num_reqs_pending > 0) {
            int req_idx = -1;
            err = MPI_Waitany(nranks-1, reqs, &req_idx, MPI_STATUS_IGNORE);
            if (err != MPI_SUCCESS) {
                op_printf("ERROR: MPI_Waitany() failed\n");
                op_exit();
                exit(EXIT_FAILURE);
            }
            int r = req_idx+1;
            if (r < 0 || r > (nranks-1)) {
                op_printf("ERROR: MPI_Waitany() has selected invalid r=%d\n", r);
                op_exit();
                exit(EXIT_FAILURE);
            }
            num_reqs_pending--;

            // printf("Root received table from rank %d:\n", r);
            // printf("----------------------------\n");
            // printf("Thread, CpuId , Event , Kernel , Level , Count\n");
            // for (int i=0; i<n_rows_each[r]; i++) {
            //     printf("R%d | %d , %d , %s , %s , %d , %lld\n", 
            //         r, 
            //         papi_tables[r][i].thread_id, 
            //         papi_tables[r][i].cpu_id, 
            //         papi_tables[r][i].event_name, papi_tables[r][i].kernel_name, 
            //         papi_tables[r][i].level, papi_tables[r][i].count);
            // }
            // printf("R%d |----------------------------\n", r);
        }
    } else {
        // printf("Rank %d sending table:\n", rank);
        // printf("----------------------------\n");
        // printf("R%d | Thread, CpuId , Event , Kernel , Level , Count\n", rank);
        // for (int i=0; i<n_rows; i++) {
        //     printf("R%d | %d , %d , %s , %s , %d , %lld\n", 
        //         rank, 
        //         my_table[i].thread_id, 
        //         my_table[i].cpu_id, 
        //         my_table[i].event_name, my_table[i].kernel_name, 
        //         my_table[i].level, my_table[i].count);
        // }
        // printf("R%d |----------------------------\n", rank);
        MPI_Request req;
        err = MPI_Isend(my_table, n_rows, MPI_PapiRowType, 0, 0, MPI_COMM_WORLD, &req);
        if (err != MPI_SUCCESS) {
            op_printf("ERROR: MPI_Isend() failed\n");
            op_exit();
            exit(EXIT_FAILURE);
        }
        MPI_Wait(&req, MPI_STATUS_IGNORE);
    }

    // Write out data:
    if (rank == 0) {
        std::string filepath = std::string(output_file_prefix);
        if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
            filepath += ".";
        }
        filepath += "PAPI.csv";

        std::ofstream outfile;
        outfile.open(filepath.c_str(), std::ios_base::out);
        std::ostringstream header;
        header << "Rank";
        header << ",Thread";
        header << ",CpuId";
        header << ",Partitioner";
        header << ",PAPI event";
        header << ",kernel";
        header << ",level";
        header << ",count";
        outfile << header.str() << std::endl;

        for (int row_idx=0; row_idx<n_rows; row_idx++) {
            std::ostringstream event_data_line;
            event_data_line << 0;
            event_data_line << "," << my_table[row_idx].thread_id;
            event_data_line << "," << my_table[row_idx].cpu_id;
            event_data_line << "," << conf.partitioner_string;
            event_data_line << "," << my_table[row_idx].event_name;
            event_data_line << "," << my_table[row_idx].kernel_name;
            event_data_line << "," << my_table[row_idx].level;
            event_data_line << "," << my_table[row_idx].count;
            outfile << event_data_line.str() << std::endl;
        }
        for (int r=1; r<nranks; r++) {
            for (int row_idx=0; row_idx<n_rows_each[r]; row_idx++) {
                std::ostringstream event_data_line;
                event_data_line << r;
                event_data_line << "," << papi_tables[r][row_idx].thread_id;
                event_data_line << "," << papi_tables[r][row_idx].cpu_id;
                event_data_line << "," << conf.partitioner_string;
                event_data_line << "," << papi_tables[r][row_idx].event_name;
                event_data_line << "," << papi_tables[r][row_idx].kernel_name;
                event_data_line << "," << papi_tables[r][row_idx].level;
                event_data_line << "," << papi_tables[r][row_idx].count;
                outfile << event_data_line.str() << std::endl;
            }
        }

        outfile.close();
    }
}
#endif // End MPI

#endif // End PAPI

#endif // End header
