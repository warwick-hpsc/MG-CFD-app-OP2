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

#ifndef PERF_FUNCS_H
#define PERF_FUNCS_H

// OP2 op_printf
#include "op_seq.h"

#include "config.h"
#include "const.h"
#include "global.h"

typedef struct perf_data_row_t
{
    char kernel_name[TABLE_STRING_LENGTH];
    int level;
    double computeTime;
    double syncTime;
    long iters;
} perf_data_row;

#ifdef _MPI
#include <mpi.h>
inline void dump_perf_data_to_file_mpi(
    perf_data_row* my_table, int n_rows, 
    char* output_file_prefix);
#endif

inline void dump_perf_data_to_file(
    #ifdef VERIFY_OP2_TIMING
        double* flux_kernel_compute_times, 
        double* flux_kernel_sync_times, 
    #endif
    long* flux_kernel_iter_counts, 
    char* output_file_prefix)
{
    int n_rows = levels;
    perf_data_row* my_table = (perf_data_row*)malloc(levels*sizeof(perf_data_row));
    for (int l=0; l<levels; l++) {
        strncpy(my_table[l].kernel_name, "compute_flux_edge_kernel", TABLE_STRING_LENGTH);
        my_table[l].level = l;
        #ifdef VERIFY_OP2_TIMING
            my_table[l].computeTime = flux_kernel_compute_times[l];
            my_table[l].syncTime    = flux_kernel_sync_times[l];
        #endif
        my_table[l].iters = flux_kernel_iter_counts[l];
    }

    #ifdef _MPI
        dump_perf_data_to_file_mpi(my_table, n_rows, output_file_prefix);
        return;
    #endif

    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += "PerfData.csv";

    std::ostringstream header;
    header << "partitioner";
    header << ",kernel";
    header << ",level";
    #ifdef VERIFY_OP2_TIMING
        header << ",computeTime";
        header << ",syncTime";
    #endif
    header << ",iters";

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::out);
    outfile << header.str() << std::endl;
    for (int l=0; l<levels; l++) {
        std::ostringstream data_line;
        data_line << conf.partitioner_string;
        data_line << ',' << my_table[l].kernel_name;
        data_line << ',' << my_table[l].level;
        #ifdef VERIFY_OP2_TIMING
            data_line << ',' << my_table[l].computeTime;
            data_line << ',' << my_table[l].syncTime;
        #endif
        data_line << ',' << my_table[l].iters;
        outfile << data_line.str() << std::endl;
    }
    outfile.close();
}

#ifdef _MPI
#include <mpi.h>
inline void dump_perf_data_to_file_mpi(
    perf_data_row* my_table, int n_rows, 
    char* output_file_prefix)
{
    // Non-root ranks send data to root, which writes to file.

    int err;

    // Create MPI_Datatype MPI_PerfRowType
    MPI_Datatype types[5] = { MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG };
    int blocklengths[5] = { TABLE_STRING_LENGTH, 1, 1, 1, 1 };
    int type_sizes[5] = { sizeof(char), sizeof(int), sizeof(double), sizeof(double), sizeof(long) };
    MPI_Aint displacements[5];
    displacements[0] = offsetof(perf_data_row, kernel_name);
    displacements[1] = offsetof(perf_data_row, level);
    displacements[2] = offsetof(perf_data_row, computeTime);
    displacements[3] = offsetof(perf_data_row, syncTime);
    displacements[4] = offsetof(perf_data_row, iters);
    MPI_Datatype MPI_PerfRowType;
    err = MPI_Type_create_struct(
        5,
        blocklengths,
        displacements,
        types,
        &MPI_PerfRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to create custom MPI 'PerfData table' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }
    err = MPI_Type_commit(&MPI_PerfRowType);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: Failed to commit custom MPI 'PerfData row' type\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    // Send/receive tables
    int nranks, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    perf_data_row** perf_data_tables = NULL;
    int* n_rows_each = (int*)malloc(nranks * sizeof(int));
    err = MPI_Gather(&n_rows, 1, MPI_INT, n_rows_each, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err != MPI_SUCCESS) {
        op_printf("ERROR: MPI_Irecv() failed\n");
        op_exit();
        exit(EXIT_FAILURE);
    }

    if (rank == 0) {
        perf_data_tables = (perf_data_row**)malloc(nranks*sizeof(perf_data_row*));
        for (int r=0; r<nranks; r++) {
            perf_data_tables[r] = (perf_data_row*)malloc(n_rows_each[r]*sizeof(perf_data_row));
        }

        MPI_Request* reqs = (MPI_Request*)malloc((nranks-1)*sizeof(MPI_Request));

        for (int r=1; r<nranks; r++) {
            err = MPI_Irecv(perf_data_tables[r], n_rows_each[r], MPI_PerfRowType, r, 0, MPI_COMM_WORLD, &reqs[r-1]);
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
            // printf("Kernel , Level , Iters\n");
            // for (int i=0; i<n_rows_each[r]; i++) {
            //     printf("R%d |%s , %d , %ld\n", 
            //         r, 
            //         perf_data_tables[r][i].kernel_name, 
            //         perf_data_tables[r][i].level, perf_data_tables[r][i].iters);
            // }
            // printf("R%d |----------------------------\n", r);
        }
    } else {
        // printf("Rank %d sending table:\n", rank);
        // printf("----------------------------\n");
        // printf("R%d | Kernel , Level , Iters\n", rank);
        // for (int i=0; i<n_rows; i++) {
        //     printf("R%d | %s , %d , %ld\n", 
        //         rank, 
        //         my_table[i].kernel_name, 
        //         my_table[i].level, my_table[i].iters);
        // }
        // printf("R%d |----------------------------\n", rank);
        MPI_Request req;
        err = MPI_Isend(my_table, n_rows, MPI_PerfRowType, 0, 0, MPI_COMM_WORLD, &req);
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
        filepath += "PerfData.csv";

        std::ostringstream header;
        header << "Rank";
        header << ",partitioner";
        header << ",kernel";
        header << ",level";
        #ifdef VERIFY_OP2_TIMING
            header << ",computeTime";
            header << ",syncTime";
        #endif
        header << ",iters";

        std::ofstream outfile;
        outfile.open(filepath.c_str(), std::ios_base::out);
        outfile << header.str() << std::endl;

        for (int l=0; l<levels; l++) {
            std::ostringstream data_line;
            data_line << 0;
            data_line << ',' << conf.partitioner_string;
            data_line << ',' << my_table[l].kernel_name;
            data_line << ',' << my_table[l].level;
            #ifdef VERIFY_OP2_TIMING
                data_line << ',' << my_table[l].computeTime;
                data_line << ',' << my_table[l].syncTime;
            #endif
            data_line << ',' << my_table[l].iters;
            outfile << data_line.str() << std::endl;
        }

        for (int r=1; r<nranks; r++) {
            for (int l=0; l<levels; l++) {
                std::ostringstream data_line;
                data_line << r,
                data_line << ',' << conf.partitioner_string;
                data_line << ',' << perf_data_tables[r][l].kernel_name;
                data_line << ',' << perf_data_tables[r][l].level;
                #ifdef VERIFY_OP2_TIMING
                    data_line << ',' << perf_data_tables[r][l].computeTime;
                    data_line << ',' << perf_data_tables[r][l].syncTime;
                #endif
                data_line << ',' << perf_data_tables[r][l].iters;
                outfile << data_line.str() << std::endl;
            }
        }

        outfile.close();
    }
}
#endif // end MPI

inline void dump_file_io_perf_data_to_file(
    double walltime, 
    double* file_io_times,
    int number_of_file_io_writes,
    char* output_file_prefix)
{
    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    #ifdef _MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank != 0) {
            // Parallel I/O time is same for each rank, so only one 
            // needs to write out.
            return;
        }
    #endif
    filepath += "FileIoTimes.csv";

    std::ostringstream header;
    header << ",partitioner";
    header << ",level";
    header << ",writeInterval";
    header << ",numberOfWrites";
    header << ",fileIoTime";
    header << ",wallTime";

    std::ofstream outfile;

    outfile.open(filepath.c_str(), std::ios_base::out);
    outfile << header.str() << std::endl;

    // for (int l=0; l<levels; l++) {
    // Currently, only benchmark file I/O performance on 
    // base level:
    for (int l=0; l<1; l++) {
        std::ostringstream data_line;

        data_line << ',' << conf.partitioner_string;
        data_line << ',' << l;

        data_line << ',' << conf.output_flow_interval;
        data_line << ',' << number_of_file_io_writes;
        data_line << ',' << file_io_times[l];

        data_line << ',' << walltime;

        outfile << data_line.str() << std::endl;
    }

    outfile.close();
}

#endif
