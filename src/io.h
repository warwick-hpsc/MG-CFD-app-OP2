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


#ifndef IO_H
#define IO_H

#include "hdf5.h"

#include "global.h"

void read_input_dat(
    const char* file_name, 
    int* problem_size, 
    int* levels, 
    int* base_array_index,
    std::string** layers, 
    std::string** mg_connectivity)
{
    std::ifstream file(file_name);

    bool have_size=false, have_num_levels=false, have_mesh_name=false, have_mesh_filenames=false, have_mg_mapping=false;
    std::string file_line;

    if (!file.is_open()) {
        fprintf(stderr, "Error: Could not open input file '%s'\n", file_name);
        DEBUGGABLE_ABORT
    }
    else {
        while (std::getline(file, file_line))
        {
            if (file_line.c_str()[0] == '#')
                continue;

            std::istringstream str_iss(file_line);
            if (file_line.c_str()[0] == '[')
            {
                if (strcmp(file_line.c_str(), "[levels]")==0)
                {
                    if (!have_num_levels)
                    {
                        fprintf(stderr, "Error parsing %s: Need to know number of levels before parsing level filenames\n", file_name);
                        DEBUGGABLE_ABORT
                    }
                    have_mesh_filenames = true;
                    *layers = new std::string[(*levels)];
                    for (int i=0; i<(*levels); i++)
                    {
                        std::string key;
                        std::string value;
                        if (!std::getline(file, file_line))
                        {
                            fprintf(stderr, "Error parsing %s: Have reached EOF before reading all mesh filenames\n", file_name);
                            DEBUGGABLE_ABORT
                        }
                        else
                        {
                            str_iss.clear();
                            str_iss.seekg(0, std::ios::beg);
                            str_iss.str(file_line);
                            if (!std::getline(str_iss, key, '='))
                            {
                                fprintf(stderr, "Error parsing '%s': Was expecting a key-value pair following [levels]\n", file_name);
                                DEBUGGABLE_ABORT
                            }
                            else
                            {
                                if (std::getline(str_iss, value)) {
                                    key = trim(key);
                                    value = trim(value);
                                    int idx = atoi(key.c_str());
                                    (*layers)[idx] = strdup(value.c_str());
                                }
                            }
                        }
                    }
                }
                else if (strcmp(file_line.c_str(), "[mg_mapping]")==0)
                {
                    if (!have_num_levels) {
                        fprintf(stderr, "Error parsing %s: Need to know number of levels before parsing level filenames\n", file_name);
                        DEBUGGABLE_ABORT
                    }
                    have_mg_mapping = true;
                    *mg_connectivity = new std::string[(*levels)-1];
                    for (int i=0; i<(*levels)-1; i++)
                    {
                        std::string key;
                        std::string value;
                        if (!std::getline(file, file_line)) {
                            fprintf(stderr, "Error parsing %s: Have reached EOF before reading all MG filenames\n", file_name);
                            DEBUGGABLE_ABORT
                        }
                        else
                        {
                            str_iss.clear();
                            str_iss.seekg(0, std::ios::beg);
                            str_iss.str(file_line);
                            if (!std::getline(str_iss, key, '=')) {
                                fprintf(stderr, "Error parsing '%s': Was expecting a key-value pair following [mg_mapping]\n", file_name);
                                DEBUGGABLE_ABORT
                            } else {
                                if (std::getline(str_iss, value)) {
                                    key = trim(key);
                                    value = trim(value);
                                    int idx = atoi(key.c_str());
                                    (*mg_connectivity)[idx] = strdup(value.c_str());
                                }
                            }
                        }
                    }
                    have_mg_mapping = true;
                }
            }
            else
            {
                std::string key;
                std::string value;
                if (std::getline(str_iss, key, '=')) {
                    // Line contains a key-value pair
                    if (std::getline(str_iss, value))
                    {
                        key = trim(key);
                        value = trim(value);

                        if (strcmp(key.c_str(), "size")==0) {
                            *problem_size = atoi(value.c_str());
                            have_size = true;
                        }
                        else if (strcmp(key.c_str(), "num_levels")==0) {
                            *levels = atoi(value.c_str());
                            have_num_levels = true;
                        }
                        else if (strcmp(key.c_str(), "base_array_index")==0) {
                            *base_array_index = atoi(value.c_str());
                        }
                        else if (strcmp(key.c_str(), "mesh_name")==0) {
                            if (strcmp(value.c_str(), "la_cascade")==0) {
                                mesh_name = MESH_LA_CASCADE;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "rotor37")==0) {
                                mesh_name = MESH_ROTOR37;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "fvcorr")==0) {
                                mesh_name = MESH_FVCORR;
                                have_mesh_name = true;
                            }
                            else if (strcmp(value.c_str(), "m6wing")==0) {
                                mesh_name = MESH_M6WING;
                                have_mesh_name = true;
                            }
                            else {
                                fprintf(stderr, "Error parsing %s: Unknown mesh_name '%s'\n", file_name, value.c_str());
                                DEBUGGABLE_ABORT
                            }
                        }
                    }
                }
            }
        }
    }

    if (!have_size) {
        fprintf(stderr, "Error parsing '%s': size not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_num_levels) {
        fprintf(stderr, "Error parsing '%s': number of levels not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_mesh_name) {
        fprintf(stderr, "Error parsing '%s': mesh name not present\n", file_name);
        DEBUGGABLE_ABORT
    }
    if (!have_mesh_filenames) {
        fprintf(stderr, "Error parsing '%s': mesh filenames not present\n", file_name);
        DEBUGGABLE_ABORT
    }

    if (!have_mg_mapping) {
        *mg_connectivity = new std::string[(*levels)-1];
        for (int i=0; i<(*levels)-1; i++) {
            (*mg_connectivity)[i] = (char*)malloc(sizeof(char));
            (*mg_connectivity)[i][0] = '\0';
        }
    }
}

void read_mg_connectivity(const char* file_name, int** mg_connectivity, int* mgc)
{
    std::ifstream file(file_name);
    *mgc = 0;
    if(file.is_open())
    {
        file >> *mgc;
        *mg_connectivity = alloc<int>(*mgc);

        for(int i = 0; i < *mgc; i++)
        {
            file >> (*mg_connectivity)[i];
        }
    }
    else
    {
        std::cout << "could not open mg file: '" << file_name << "'" << std::endl;
        DEBUGGABLE_ABORT
    }
}

void dump_double_array(
    double* array, 
    int n, 
    int dim, 
    const char* name,
    std::string filenamePrefix, 
    std::string filenameSuffix)
{
    std::string file_name;
    if (filenamePrefix != "") file_name += filenamePrefix;
    file_name += std::string(name);
    if (filenameSuffix != "") file_name += filenameSuffix;

    op_printf("Dumping '%s' to '%s'\n", name, file_name.c_str());

    FILE *file = fopen(file_name.c_str(), "w");
    for (int i=0; i<n; i++) {
        fprintf(file, "%.17e", array[i*dim]);
        for (int j=1; j<dim; j++) {
            fprintf(file, " %.17e", array[i*dim+j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void dump_int_array(
    int* array, 
    int n, 
    int dim, 
    const char* name,
    std::string filenamePrefix, 
    std::string filenameSuffix)
{
    std::string file_name;
    if (filenamePrefix != "") file_name += filenamePrefix;
    file_name += std::string(name);
    if (filenameSuffix != "") file_name += filenameSuffix;

    FILE *file = fopen(file_name.c_str(), "w");
    for (int i=0; i<n; i++) {
        for (int j=0; j<dim; j++) {
            fprintf(file, " %d", array[i*dim+j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void dump_edge_neighbours(
    edge_neighbour* edges, 
    int nedges, 
    const char* name,
    std::string filenamePrefix, 
    std::string filenameSuffix)
{
    std::string file_name;
    if (filenamePrefix != "") file_name += filenamePrefix;
    file_name += std::string(name);
    if (filenameSuffix != "") file_name += filenameSuffix;

    FILE *file = fopen(file_name.c_str(), "w");
    for (int i=0; i<nedges; i++) {
        fprintf(file, "%d %d %.17e %.17e %.17e\n", edges[i].a, edges[i].b, edges[i].x, edges[i].y, edges[i].z);
    }
    fclose(file);
}

inline void dump_perf_data_to_file(
    int rank, 
    int num_levels, 
    #ifdef VERIFY_OP2_TIMING
        double* flux_kernel_compute_times, 
        double* flux_kernel_sync_times, 
    #endif
    long* flux_kernel_iter_counts, 
    char* output_file_prefix)
{
    std::string filepath = std::string(output_file_prefix);
    if (filepath.length() > 1 && filepath.at(filepath.size()-1) != '/') {
        filepath += ".";
    }
    filepath += std::string("P=") + number_to_string(rank);
    filepath += ".PerfData.csv";

    bool write_header = false;
    std::ifstream f(filepath.c_str());
    if (!f || f.peek() == std::ifstream::traits_type::eof()) {
        write_header = true;
    }
    f.close();

    std::ostringstream header;

    if (write_header) {
        header << "rank";
        header << ",partitioner";
        header << ",kernel";
        header << ",level";
        #ifdef VERIFY_OP2_TIMING
            header << ",computeTime";
            header << ",syncTime";
        #endif
        header << ",iters";
    }

    std::ofstream outfile;
    outfile.open(filepath.c_str(), std::ios_base::app);
    if (write_header) outfile << header.str() << std::endl;

    for (int l=0; l<num_levels; l++) {
        std::ostringstream data_line;

        data_line << rank;
        data_line << ',' << conf.partitioner_string;
        data_line << ",compute_flux_edge_kernel";
        data_line << ',' << l;

        #ifdef VERIFY_OP2_TIMING
            data_line << ',' << flux_kernel_compute_times[l];
            data_line << ',' << flux_kernel_sync_times[l];
        #endif
        data_line << ',' << flux_kernel_iter_counts[l];

        outfile << data_line.str() << std::endl;
    }

    outfile.close();
}

#endif

