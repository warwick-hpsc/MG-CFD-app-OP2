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
#include <regex>

// OP2:
#include "op_seq.h"

#include "config.h"
#include "const.h"
#include "global.h"

extern int n_events;
extern char** event_names;
extern int likwid_gid;
extern CpuTopology_t likwid_topo;

extern int n_cpus_monitored;
extern int last_cpu_id;

#define long_long long long
extern long_long** flux_kernel_event_counts;
extern long_long** ustream_kernel_event_counts;

#define TABLE_STRING_LENGTH 128
typedef struct likwid_row_t
{
    int thread_id;
    int cpu_id;
    char event_name[TABLE_STRING_LENGTH];
    char kernel_name[TABLE_STRING_LENGTH];
    int level;
    long long int count;
} likwid_row;

#ifdef _MPI
#include <mpi.h>
inline void dump_likwid_counters_to_file_mpi(
    int rank, 
    likwid_row* my_table, int n_rows, 
    char* output_file_prefix);
#endif

inline int likwid_get_num_threads() {
    #ifdef _OMP
        return omp_get_max_threads();
    #else
        return 1;
    #endif
}

inline void clear_likwid()
{
    perfmon_finalize();
    affinity_finalize();
    topology_finalize();
}

inline void my_likwid_start()
{
    if (n_events == 0) {
        return;
    }

    if (last_cpu_id == -1) {
        last_cpu_id = sched_getcpu();
    } else {
        if (last_cpu_id != sched_getcpu()) {
            op_printf("ERROR: Process has moved to a different processor (%d -> %d), this invalidates Likwid counts. Try again with pinning.\n", sched_getcpu(), last_cpu_id);
            clear_likwid();
            op_exit();
            exit(EXIT_FAILURE);
        }
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
		// for (int t=0; t<likwid_topo->numHWThreads; t++) {
        // for (int t=0; t<likwid_get_num_threads(); t++) {
        for (int t=0; t<n_cpus_monitored; t++) {
			double result_f = perfmon_getResult(likwid_gid, e, t);
            // printf("perfmon_getResult(gid=%d, e=%d, t=%d) = %.2e (current_level=%d)\n", likwid_gid, e, t, result_f, current_level);
			// long_long result = (long_long)result_f;
			// event_counts[t][current_level*n_events + e] += result;
            event_counts[t][current_level*n_events + e] += result_f;
            // printf("event_counts[t][current_level*n_events + e] = %f\n", event_counts[t][current_level*n_events + e]);
            // event_counts[t][current_level*n_events + e] = result_f;
		}
	}
}

inline void init_likwid()
{
    int err, j;
    err = topology_init();
    if (err < 0) {
        op_printf("Failed to initialize LIKWID's topology module\n");
        op_exit();
        exit(EXIT_FAILURE);
    }
    CpuInfo_t info = get_cpuInfo();
    likwid_topo = get_cpuTopology();
    affinity_init();

    // n_cpus_monitored = likwid_topo->numHWThreads;
    // int* cpus = (int*)malloc(n_cpus_monitored * sizeof(int));
    // for (j=0; j<n_cpus_monitored;j++) {
    //     cpus[j] = likwid_topo->threadPool[j].apicId;
    // }

    #ifdef _OMP
        n_cpus_monitored = likwid_get_num_threads();
        int* cpus = (int*)malloc(n_cpus_monitored * sizeof(int));
        #pragma omp parallel
        {
            cpus[omp_get_thread_num()] = sched_getcpu();
        }
        // As we know thread placement, can check that none
        // are on same CPU
        for (int t1=0; t1<(likwid_get_num_threads()-1); t1++) {
            for (int t2=t1+1; t2<likwid_get_num_threads(); t2++) {
                if (cpus[t1] == cpus[t2]) {
                    op_printf("ERROR: Threads %d and %d are on same CPU %d\n", t1, t2, cpus[t1]);
                    op_exit();
                    exit(EXIT_FAILURE);
                }
            }
        }
    #else
        n_cpus_monitored = 1;
        int* cpus = (int*)malloc(n_cpus_monitored * sizeof(int));
        cpus[0] = sched_getcpu();
    #endif

    err = perfmon_init(n_cpus_monitored, cpus);
    if (err < 0) {
        op_printf("Failed to initialize LIKWID's performance monitoring module\n");
        op_exit();
        topology_finalize();
        exit(EXIT_FAILURE);
    }

    // numa_init();
    // perfmon_init_maps();
    // perfmon_check_counter_map(0);

    n_events = 0;
    likwid_gid = 0;
    last_cpu_id = -1;
}

inline bool uncore_monitoring_permitted()
{
    // Likwid library does not prevent processes in same socket read/writing same
    // uncore MSRs. So need to prevent this manually, somewhat like likwid-mpirun does:

    #ifndef _MPI
        // No intra-process contention for uncore, and Likwid library
        // takes care of intra-thread contention, so safe.
        return true;
    #else
        bool permit_uncore_access = true;
        int err;

        //
        // Determine if this rank has lowest-count CPU among ranks in that socket:
        //

        // Create communicator for ranks in same node
        MPI_Comm node_comm;
        err =  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, 0, &node_comm);
        if (err != MPI_SUCCESS) {
            op_printf("ERROR: Failed to create per-node communicator\n");
            op_exit();
            exit(EXIT_FAILURE);
        }

        int node_size = 0;
        MPI_Comm_size(node_comm, &node_size);
        if (node_size == 0) {
            op_printf("ERROR: MPI_Comm_size() has returned 0 for my node comm\n");
            op_exit();
            exit(EXIT_FAILURE);
        }

        int my_node_rank = -1;
        MPI_Comm_rank(node_comm, &my_node_rank);
        if (my_node_rank == -1) {
            op_printf("ERROR: MPI_Comm_rank() has returned -1 for my node rank\n");
            op_exit();
            exit(EXIT_FAILURE);
        }

        // Get occupied cpus
        uint my_cpu = sched_getcpu();
        uint* node_occupied_cpus = (uint*)malloc(node_size*sizeof(uint));
        for (int i=0; i<node_size; i++) node_occupied_cpus[i] = 999;
        node_occupied_cpus[my_node_rank] = sched_getcpu();
        err = MPI_Allgather(&my_cpu, 1, MPI_INT, node_occupied_cpus, 1, MPI_INT, node_comm);

        // Map cpus to sockets
        uint num_cpus = likwid_topo->numCoresPerSocket * likwid_topo->numSockets * likwid_topo->numThreadsPerCore;
        uint* cpu_to_skt = (uint*)malloc(num_cpus * sizeof(uint));
        for (uint i=0; i<num_cpus; i++) cpu_to_skt[i] = 99;
        for (uint j=0; j<likwid_topo->numHWThreads; j++) {
            uint apicId = likwid_topo->threadPool[j].apicId;
            uint threadId = likwid_topo->threadPool[j].threadId;
            uint coreId = likwid_topo->threadPool[j].coreId;
            uint packageId = likwid_topo->threadPool[j].packageId;
            uint cpuId = apicId;
            if (cpu_to_skt[cpuId] == 99) {
                cpu_to_skt[cpuId] = packageId;
            }
            else {
                if (cpu_to_skt[cpuId] != packageId) {
                    op_printf("ERROR: Likwid CpuTopology is corrupt, mapping CPU %d to different sockets (%d, %d)\n", 
                        cpuId, 
                        cpu_to_skt[cpuId], 
                        packageId);
                    op_exit();
                    clear_likwid();
                    exit(EXIT_FAILURE);
                }
            }
        }

        uint my_skt = cpu_to_skt[my_cpu];
        for (int i=0; i<node_size; i++) {
            const uint cpu = node_occupied_cpus[i];
            if (cpu_to_skt[cpu] == my_skt && cpu < my_cpu) {
                permit_uncore_access = false;
            }
        }

        return permit_uncore_access;
    #endif
}

inline void load_likwid_events()
{
    int err;
    #ifdef _MPI

    #endif

    bool permit_uncore_access = uncore_monitoring_permitted();

    std::vector<std::string> events;
    std::regex likwid_event_rgx("([A-Z0-9_]*):([A-Z0-9_]*)");
    std::regex likwid_mbox_rgx("MBOX.*");

    std::string line;
    std::ifstream file_if(conf.papi_config_file);
    if(!file_if.is_open()) {
        op_printf("WARNING: Failed to open Likwid config: '%s'\n", conf.papi_config_file);
        op_exit();
        exit(EXIT_FAILURE);
    } else {
        while(std::getline(file_if, line)) {
            if (line.c_str()[0] == '#' || strcmp(line.c_str(), "")==0) {
                continue;
            }

            std::smatch match;
            if (!std::regex_match(line, match, likwid_event_rgx)) {
                op_printf("ERROR: Provided Likwid event is malformed: %s\n", line.c_str());
                op_exit();
                clear_likwid();
                exit(EXIT_FAILURE);
            }
            std::string rhs = match.str(2);

            if (std::regex_match(rhs, likwid_mbox_rgx)) {
                if (permit_uncore_access) {
                    events.push_back(line);
                }
            } else {
                events.push_back(line);
            }
        }
    }

    n_events = events.size();
    if (n_events > 0) {
        int j, n;

        const int nt = n_cpus_monitored;
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
            n_events = 0;
            return;
        }
        err = perfmon_setupCounters(likwid_gid);
        if (err < 0) {
            op_printf("Failed to setup group %d in LIKWID's performance monitoring module\n", likwid_gid);
            n_events = 0;
            return;
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
    const int nt = n_cpus_monitored;

    int n_rows = nt * n_events * levels;
    if (conf.measure_mem_bound) {
        n_rows *= 2;
    }
    likwid_row* my_table = (likwid_row*)malloc(n_rows*sizeof(likwid_row));
    int* cpu_ids = (int*)malloc(n_cpus_monitored*sizeof(int));
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
        for (int eid=0; eid<n_events; eid++) {
            for (int l=0; l<levels; l++) {
                my_table[row_idx].thread_id = tid;
                my_table[row_idx].cpu_id = cpu_ids[tid];
                strncpy(my_table[row_idx].event_name, event_names[eid], TABLE_STRING_LENGTH);
                strncpy(my_table[row_idx].kernel_name, "compute_flux_edge_kernel", TABLE_STRING_LENGTH);
                my_table[row_idx].level = l;
                const int idx = l*n_events + eid;
                my_table[row_idx].count = flux_kernel_event_counts[tid][idx];
                row_idx++;
            }

            if (conf.measure_mem_bound) {
                for (int l=0; l<levels; l++) {
                    my_table[row_idx].thread_id = tid;
                    my_table[row_idx].cpu_id = cpu_ids[tid];
                    strncpy(my_table[row_idx].event_name, event_names[eid], TABLE_STRING_LENGTH);
                    strncpy(my_table[row_idx].kernel_name, "unstructured_stream_kernel", TABLE_STRING_LENGTH);
                    my_table[row_idx].level = l;
                    const int idx = l*n_events + eid;
                    my_table[row_idx].count = ustream_kernel_event_counts[tid][idx];
                    row_idx++;
                }
            }
        }
    }

    #ifdef _MPI
        dump_likwid_counters_to_file_mpi(rank, my_table, n_rows, output_file_prefix);
        return;
    #endif

    // Write out data:
    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "Likwid.csv";

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::out);
    std::ostringstream header;
    header << "Rank";
    header << ",Thread";
    header << ",CpuId";
    header << ",Partitioner";
    header << ",Likwid event";
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

    outfile.close();
}

// Functionality below should be moved into OP2
#ifdef _MPI
#include <mpi.h>
inline void dump_likwid_counters_to_file_mpi(
	int rank, 
    likwid_row* my_table, int n_rows, 
    char* output_file_prefix)
{
    // Non-root ranks send data to root, which writes to file.

    int err;

    // Create MPI_Datatype MPI_LikwidRowType
    int types[6] = { MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_LONG_LONG_INT};
    int blocklengths[6] = { 1, 1, TABLE_STRING_LENGTH, TABLE_STRING_LENGTH, 1, 1 };
    int type_sizes[6] = { sizeof(int), sizeof(int), sizeof(char), sizeof(char), sizeof(int), sizeof(long long int) };
    MPI_Aint displacements[6];
    displacements[0] = offsetof(likwid_row, thread_id);
    displacements[1] = offsetof(likwid_row, cpu_id);
    displacements[2] = offsetof(likwid_row, event_name);
    displacements[3] = offsetof(likwid_row, kernel_name);
    displacements[4] = offsetof(likwid_row, level);
    displacements[5] = offsetof(likwid_row, count);
    MPI_Datatype MPI_LikwidRowType;
    err = MPI_Type_create_struct(
        6,
        blocklengths,
        displacements,
        types,
        &MPI_LikwidRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to create custom MPI 'Likwid table' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }
    err = MPI_Type_commit(&MPI_LikwidRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to commit custom MPI 'Likwid row' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    // Send/receive tables
    int nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    likwid_row** likwid_tables = NULL;
    int* n_rows_each = (int*)malloc(nranks * sizeof(int));
    err = MPI_Gather(&n_rows, 1, MPI_INT, n_rows_each, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: MPI_Irecv() failed\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        likwid_tables = (likwid_row**)malloc(nranks*sizeof(likwid_row*));
        for (int r=0; r<nranks; r++) {
            likwid_tables[r] = (likwid_row*)malloc(n_rows_each[r]*sizeof(likwid_row));
        }

        MPI_Request* reqs = (MPI_Request*)malloc((nranks-1)*sizeof(MPI_Request));

        for (int r=1; r<nranks; r++) {
            err = MPI_Irecv(likwid_tables[r], n_rows_each[r], MPI_LikwidRowType, r, 0, MPI_COMM_WORLD, &reqs[r-1]);
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
            //         likwid_tables[r][i].thread_id, 
            //         likwid_tables[r][i].cpu_id, 
            //         likwid_tables[r][i].event_name, likwid_tables[r][i].kernel_name, 
            //         likwid_tables[r][i].level, likwid_tables[r][i].count);
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
        err = MPI_Isend(my_table, n_rows, MPI_LikwidRowType, 0, 0, MPI_COMM_WORLD, &req);
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
        filepath += "Likwid.csv";

        std::ofstream outfile;
        outfile.open(filepath.c_str(), std::ios_base::out);
        std::ostringstream header;
        header << "Rank";
        header << ",Thread";
        header << ",CpuId";
        header << ",Partitioner";
        header << ",Likwid event";
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
                event_data_line << "," << likwid_tables[r][row_idx].thread_id;
                event_data_line << "," << likwid_tables[r][row_idx].cpu_id;
                event_data_line << "," << conf.partitioner_string;
                event_data_line << "," << likwid_tables[r][row_idx].event_name;
                event_data_line << "," << likwid_tables[r][row_idx].kernel_name;
                event_data_line << "," << likwid_tables[r][row_idx].level;
                event_data_line << "," << likwid_tables[r][row_idx].count;
                outfile << event_data_line.str() << std::endl;
            }
        }

        outfile.close();
    }

    op_printf("dump_likwid_counters_to_file() complete\n");
}
#endif

#endif
