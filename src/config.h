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
    { "config_filepath",    required_argument, NULL, 'c'}, 
    { "legacy_mode",        no_argument,       NULL, 'l' },
    { "input-file",         required_argument, NULL, 'i' },
    { "input-directory",    required_argument, NULL, 'd' },
    { "papi_config_file",   required_argument, NULL, 'p' },
    { "output-file-prefix", required_argument, NULL, 'o' },
    { "num-cycles",         required_argument, NULL, 'g' },
    { "partitioner",        required_argument, NULL, 'm' },
    { "validate",           no_argument,       NULL, 'v' },
    { "output-variables",   no_argument,       (int*)&conf.output_variables,    1 },
};
#define GETOPTS "c:li:d:p:o:g:m:v"

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

    conf.num_cycles = 25;

    conf.partitioner = Partitioners::Parmetis;

    conf.output_volumes = false;
    conf.output_step_factors = false;
    conf.output_edge_mx = false;
    conf.output_edge_my = false;
    conf.output_edge_mz = false;
    conf.output_edge_p  = false;
    conf.output_edge_pe = false;
    conf.output_fluxes  = false;
    conf.output_variables = false;
    conf.output_anything = false;
}

inline void set_config_param(const char* const key, const char* const value) {
    if (strcmp(key,"config_filepath")==0) {
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
    else if (strcmp(key,"papi_config_file")==0) {
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
        conf.partitioner_string = strdup(value);
    }

    else if (strcmp(key,"output_volumes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_volumes = true;
        }
    }
    else if (strcmp(key,"output_step_factors")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_step_factors = true;
        }
    }
    else if (strcmp(key,"output_edge_fluxes")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_mx = true;
            conf.output_edge_my = true;
            conf.output_edge_mz = true;
            conf.output_edge_p  = true;
            conf.output_edge_pe = true;
        }
    }
    else if (strcmp(key,"output_edge_mx")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_mx = true;
        }
    }
    else if (strcmp(key,"output_edge_my")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_my = true;
        }
    }
    else if (strcmp(key,"output_edge_mz")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_mz = true;
        }
    }
    else if (strcmp(key,"output_edge_p")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_p = true;
        }
    }
    else if (strcmp(key,"output_edge_pe")==0) {
        if (strcmp(value, "Y")==0) {
            conf.output_edge_pe = true;
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
        set_config_param("papi_config_file", (std::string(conf.input_file_directory) + "/" + conf.papi_config_file).c_str());
    }
    #endif
}

inline bool parse_arguments(int argc, char** argv) {
    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch(optc) {
            // case 'h':
            //     print_help();
            //     return false;
            //     break;
            case 'i':
                set_config_param("input_file", strdup(optarg));
                break;
            case 'c':
                set_config_param("config_filepath", strdup(optarg));
                read_config();
                break;
            case 'd':
                set_config_param("input_file_directory", strdup(optarg));
                break;
            case 'p':
                set_config_param("papi_config_file", strdup(optarg));
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

    return true;
}

#endif
