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

#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
#include <unistd.h>

#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

#include "utils.h"

namespace Partitioners
{
    enum Partitioners {
        Inertial, 
        Parmetis, 
        Ptscotch
    };
}

namespace PartitionerMethods
{
    enum PartitionerMethods {
        Geom, 
        KWay, 
        GeomKWay,
        NotSet
    };
}


typedef struct {
    char* config_filepath;
    char* input_file;
    char* input_file_directory;
    char* output_file_prefix;

    #ifdef PAPI
    char* papi_config_file;
    #endif

    bool legacy_mode;

    int num_cycles;

    Partitioners::Partitioners partitioner;
    char* partitioner_string;

    PartitionerMethods::PartitionerMethods partitioner_method;
    char* partitioner_method_string;

    bool validate_result;

    bool output_volumes;
    bool output_step_factors;
    bool output_edge_mx;
    bool output_edge_my;
    bool output_edge_mz;
    bool output_edge_p;
    bool output_edge_pe;
    bool output_fluxes;
    bool output_variables;

    bool output_anything;
} config;

extern config conf;

static struct option long_opts[] = 
{
    { "help",               no_argument,       NULL, 'h' },
    { "config-filepath",    required_argument, NULL, 'c'}, 
    { "legacy-mode",        no_argument,       NULL, 'l' },
    { "input-file",         required_argument, NULL, 'i' },
    { "input-directory",    required_argument, NULL, 'd' },
    { "papi-config-file",   required_argument, NULL, 'p' },
    { "output-file-prefix", required_argument, NULL, 'o' },
    { "num-cycles",         required_argument, NULL, 'g' },
    { "partitioner",        required_argument, NULL, 'm' },
    { "partitioner-method", required_argument, NULL, 'r' },
    { "validate",           no_argument,       NULL, 'v' },
    { "output-variables",   no_argument,       (int*)&conf.output_variables,    1 },
    { "output-fluxes",      no_argument,       (int*)&conf.output_fluxes,       1 },
    { "output-step-factors",no_argument,       (int*)&conf.output_step_factors, 1 },
};
#define GETOPTS "hc:li:d:p:o:g:m:r:v"

inline void set_config_defaults() {
    conf.config_filepath = (char*)malloc(sizeof(char));
    conf.config_filepath[0] = '\0';
    conf.input_file = (char*)malloc(sizeof(char));
    conf.input_file[0] = '\0';
    conf.input_file_directory = (char*)malloc(sizeof(char));
    conf.input_file_directory[0] = '\0';
    conf.output_file_prefix = (char*)malloc(sizeof(char));
    conf.output_file_prefix[0] = '\0';

    #ifdef PAPI
    conf.papi_config_file = (char*)malloc(sizeof(char));
    conf.papi_config_file[0] = '\0';
    #endif

    conf.legacy_mode = false;

    conf.validate_result = false;

    conf.num_cycles = 10;

    conf.partitioner = Partitioners::Parmetis;
    conf.partitioner_method = PartitionerMethods::NotSet;
    conf.partitioner_string = (char*)malloc(sizeof(char));
    conf.partitioner_string[0] = '\0';

    conf.output_step_factors = false;
    conf.output_fluxes  = false;
    conf.output_variables = false;
    conf.output_anything = false;
}

inline void set_config_param(const char* const key, const char* const value) {
    if (strcmp(key,"config-filepath")==0) {
        conf.config_filepath = strdup(value);
    }
    else if (strcmp(key,"input_file")==0) {
        conf.input_file = strdup(value);
    }
    else if (strcmp(key,"input_file_directory")==0) {
        conf.input_file_directory = strdup(value);
    }
    else if (strcmp(key,"output_file_prefix")==0) {
        conf.output_file_prefix = strdup(value);
    }

    #ifdef PAPI
    else if (strcmp(key,"papi-config-file")==0) {
        conf.papi_config_file = strdup(value);
    }
    #endif

    else if (strcmp(key,"validate_result")==0) {
        if (strcmp(value, "Y")==0) {
            conf.validate_result = true;
        }
    }

    else if (strcmp(key, "cycles")==0) {
        conf.num_cycles = atoi(value);
    }

    else if (strcmp(key, "partitioner")==0) {
        if (strcmp(value, "inertial")==0) {
            conf.partitioner = Partitioners::Inertial;
        }
        else if (strcmp(value, "parmetis")==0) {
            conf.partitioner = Partitioners::Parmetis;
        }
        else if (strcmp(value, "ptscotch")==0) {
            conf.partitioner = Partitioners::Ptscotch;
        }
        else {
            printf("WARNING: Unknown value '%s' encountered for key '%s' during parsing of config file.\n", value, key);
        }
    }

    else if (strcmp(key, "partitioner-method")==0) {
        if (strcmp(value, "geom")==0) {
            conf.partitioner_method = PartitionerMethods::Geom;
        }
        else if (strcmp(value, "kway")==0) {
            conf.partitioner_method = PartitionerMethods::KWay;
        }
        else if (strcmp(value, "geomkway")==0) {
            conf.partitioner_method = PartitionerMethods::GeomKWay;
        }
        else {
            printf("WARNING: Unknown value '%s' encountered for key '%s' during parsing of config file.\n", value, key);
        }
    }

    else if (strcmp(key,"output_step_factors")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_step_factors = true;
        }
    }
    else if (strcmp(key,"output_fluxes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_fluxes = true;
        }
    }
    else if (strcmp(key,"output_variables")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_variables = true;
        }
    }
    else {
        printf("WARNING: Unknown key '%s' encountered during parsing of config file.\n", key);
    }
}

inline void read_config() {
    if (access(conf.config_filepath, F_OK) == -1) {
        fprintf(stderr, "ERROR: \"%s\" does not exist.\n", conf.config_filepath); fflush(stderr);
        exit(-1);
    }

    // set_config_defaults();

    std::ifstream file(conf.config_filepath);
    std::string str;
    while(std::getline(file, str))
    {
        if (str.c_str()[0] == '#')
            continue;

        #ifdef LOG
            printf("Processing line: %s\n", str.c_str());
            fflush(stdout);
        #endif

        std::istringstream str_iss(str);
        std::string key;
        if (std::getline(str_iss, key, '=')) {
            std::string value;
            if (std::getline(str_iss, value)) {
                key=trim(key);
                value=trim(value);
                set_config_param(key.c_str(), value.c_str());
            }
        }
    }

    std::string config_dirpath;
    const size_t last_slash_idx = std::string(conf.config_filepath).rfind('/');
    // printf("last_slash_idx: %d\n", last_slash_idx);
    if (last_slash_idx != std::string::npos) {
        config_dirpath = std::string(conf.config_filepath).substr(0, last_slash_idx);
    } else {
        config_dirpath = std::string("");
    }
    // printf("config_dirpath:\n");
    // printf("  '%s'\n", config_dirpath.c_str());

    // if (conf.input_file_directory[0] != '/' && config_dirpath != std::string("")) {
    //     // 'input_file_directory' is currently relative to config, need to prepend config file location:
    //     std::string input_file_directory_corrected = std::string(config_dirpath) + "/" + conf.input_file_directory;
    //     // printf("input_file_directory_corrected:\n");
    //     // printf("  %s\n", input_file_directory_corrected.c_str());
    //     set_config_param(&conf, "input_file_directory", input_file_directory_corrected.c_str());
    // }

    #ifdef PAPI
    if (conf.papi_config_file != std::string("") && conf.papi_config_file[0] != '/') {
        set_config_param("papi-config-file", (std::string(conf.input_file_directory) + "/" + conf.papi_config_file).c_str());
    }
    #endif
}

inline void print_help(void)
{
    fprintf(stderr, "MG-CFD OP2 instructions\n\n");
    fprintf(stderr, "Usage: mgcfd_* [OPTIONS] -i input.dat \n");
    fprintf(stderr, "       mgcfd_* [OPTIONS] -c conf.conf \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-h, --help    Print help\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CRITICAL ARGUMENTS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  One of these must be set:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-i, --input-file=FILEPATH\n");
    fprintf(stderr, "        multigrid input grid (.dat file)\n");
    fprintf(stderr, "-c, --config-filepath=FILEPATH\n");
    fprintf(stderr, "        config file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-d, --input-directory=DIRPATH\n");
    fprintf(stderr, "        directory path to input files\n");
    fprintf(stderr, "-o, --output-file-prefix=STRING\n");
    fprintf(stderr, "        string to prepend to output filenames\n");
    fprintf(stderr, "-p, --papi-config-file=FILEPATH\n");
    fprintf(stderr, "        file containing list of PAPI events to monitor\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-g, --num-cycles=INT\n");
    fprintf(stderr, "        number of multigrid V-cycles to perform\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-m, --partitioner=STRING\n");
    fprintf(stderr, "        specify which partitioner to use:\n");
    fprintf(stderr, "          parmetis (default)\n");
    fprintf(stderr, "          ptscotch\n");
    fprintf(stderr, "          inertial\n");
    fprintf(stderr, "-r, --partitioner-method=STRING\n");
    fprintf(stderr, "        specify which partitioner method to use. Currently only\n");
    fprintf(stderr, "        applicable to ParMETIS.\n");
    fprintf(stderr, "          geom (default)\n");
    fprintf(stderr, "          kway\n");
    fprintf(stderr, "          geomkway\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "-v, --validate-result\n");
    fprintf(stderr, "        check final state against pre-calculated solution\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "DEBUGGING ARGUMENTS\n");
    fprintf(stderr, "--output-variables\n");
    fprintf(stderr, "        write Euler equation variable values to HDF5 file\n");
    fprintf(stderr, "--output-fluxes\n");
    fprintf(stderr, "        write flux accumulations to HDF5 file\n");
    fprintf(stderr, "--output-step-factors\n");
    fprintf(stderr, "        write time-step factors to HDF5 file\n");
    fprintf(stderr, "\n");
}

inline bool parse_arguments(int argc, char** argv) {
    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch(optc) {
            case 'h':
                print_help();
                return false;
            case 'i':
                set_config_param("input_file", strdup(optarg));
                break;
            case 'c':
                set_config_param("config-filepath", strdup(optarg));
                read_config();
                break;
            case 'd':
                set_config_param("input_file_directory", strdup(optarg));
                break;
            case 'p':
                set_config_param("papi-config-file", strdup(optarg));
                break;
            case 'o':
                set_config_param("output_file_prefix", strdup(optarg));
                break;
            case 'g':
                conf.num_cycles = atoi(optarg);
                break;
            case 'm':
                set_config_param("partitioner", strdup(optarg));
                break;
            case 'r':
                set_config_param("partitioner-method", strdup(optarg));
                break;
            case 'l':
                conf.legacy_mode = true;
                break;
            case 'v':
                conf.validate_result = true;
                break;
            case '\0':
                break;
            default:
                printf("Unknown command line parameter '%c'\n", optc);
        }
    }

    conf.output_anything = 
        conf.output_volumes | conf.output_step_factors | conf.output_edge_mx | 
        conf.output_edge_my | conf.output_edge_mz | conf.output_edge_p  | 
        conf.output_edge_pe | conf.output_fluxes  | conf.output_variables;

    if (conf.partitioner_method == PartitionerMethods::NotSet) {
        // Set to method partitioner-specific default:
        if (conf.partitioner == Partitioners::Ptscotch) {
            conf.partitioner_method = PartitionerMethods::KWay;
        }
        else if (conf.partitioner == Partitioners::Parmetis) {
            conf.partitioner_method = PartitionerMethods::Geom;
        }
        else if (conf.partitioner == Partitioners::Inertial) {
            conf.partitioner_method = PartitionerMethods::Geom;
        }
    }

    // Ensure partitioner and method are compatible:
    if (conf.partitioner == Partitioners::Ptscotch) {
        if (conf.partitioner_method != PartitionerMethods::KWay) {
            printf("Incompatible method requested for PT-SCOTCH, overriding with KWay\n");
            conf.partitioner_method = PartitionerMethods::KWay;
        }
    }
    if (conf.partitioner == Partitioners::Inertial) {
        if (conf.partitioner_method != PartitionerMethods::Geom) {
            printf("Incompatible method requested for Inertial, overriding with Geom\n");
            conf.partitioner_method = PartitionerMethods::Geom;
            // Although OP2 ignores the partitioner-method argument, having MG-CFD 
            // accept them would imply to the user that it has an effect. Inertial is 
            // a geometric partitioner, so pretend it requires Geom method.
        }
    }

    // Generate partitioner string:
    int new_str_len = 0;
    if (conf.partitioner == Partitioners::Inertial || 
        conf.partitioner == Partitioners::Parmetis || 
        conf.partitioner == Partitioners::Ptscotch) {
        new_str_len += 8;
    }
    new_str_len += 1;
    if (conf.partitioner_method == PartitionerMethods::Geom || 
        conf.partitioner_method == PartitionerMethods::KWay) {
        new_str_len += 4;
    }
    else if (conf.partitioner_method == PartitionerMethods::GeomKWay) {
        new_str_len += 8;
    }
    free(conf.partitioner_string);
    conf.partitioner_string = (char*)malloc(sizeof(char)*new_str_len+1);
    int i = 0;
    switch (conf.partitioner) {
        case Partitioners::Inertial:
            strcpy(conf.partitioner_string, "inertial");
            i += 8;
            break;
        case Partitioners::Parmetis:
            strcpy(conf.partitioner_string, "parmetis");
            i += 8;
            break;
        case Partitioners::Ptscotch:
            strcpy(conf.partitioner_string, "ptscotch");
            i += 8;
            break;
    }
    strcpy(conf.partitioner_string+i, "-");
    i += 1;
    switch (conf.partitioner_method) {
        case PartitionerMethods::Geom:
            strcpy(conf.partitioner_string+i, "geom");
            i += 4;
            break;
        case PartitionerMethods::KWay:
            strcpy(conf.partitioner_string+i, "kway");
            i += 4;
            break;
        case PartitionerMethods::GeomKWay:
            strcpy(conf.partitioner_string+i, "geomkway");
            i += 8;
            break;
        case PartitionerMethods::NotSet:
            strcpy(conf.partitioner_string+i, "NotSet");
            i += 8;
            break;
    }
    strcpy(conf.partitioner_string+i, "\0");

    return true;
}

#endif

