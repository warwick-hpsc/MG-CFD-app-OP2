#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <chrono>
#include <algorithm>
#include "../src_op/const_op.h"
#include "../src/structures.h"
#include "../src_op/coupler_config.h"

void read_inputs(int *indx, int *indy, int *npx, int *npy){
	FILE *pic_inputs = fopen("pbpic.input", "r");
	int buff_len = 1024;
	char line_buff[buff_len];
	long long num_particles;
	long long num_cells;
	
	if(pic_inputs == NULL){
		fprintf(stderr, "Can't open input file pbpic.input\n");
        exit(1);
    }
	while(fgets(line_buff, buff_len, pic_inputs)){
		char* token = strtok(line_buff, " ");
		if(strcmp(token, "particles") == 0){
			num_particles = atoll(strtok(NULL, " "));
		}else if(strcmp(token, "cells") == 0){
            num_cells = atoll(strtok(NULL, " "));
		}
	}
	num_particles *= 4;
    num_cells *= 8;
	int x_cells = round(sqrt(num_cells));
	*indx = round(log(std::max(x_cells*0.01,512.0))/log(2));
	//*indx = round(log(x_cells)/log(2));
	int guess_x = pow(2,*indx);
	int y_cells = round(num_cells/guess_x);
	*indy = round(log(y_cells)/log(2));
	*npx = round(std::max(sqrt(num_particles)*0.01,256.0));
	//*npx = round(sqrt(num_particles));
	*npy = round(num_particles / *npx);
}


int find_unit(struct unit units[], struct locators relative_positions[],
				   int world_rank, char type){
	int our_unit_num = relative_positions[world_rank].placelocator;
    int unit_count = 0;
    int work_count = 1; //since units start from 1
    bool found = false;
	if( type == 'C'){
		while(!found){
            if(units[unit_count].type == 'C' && our_unit_num == work_count){
                found=true;
            }else{
                if(units[unit_count].type == 'C'){
                    work_count++;
                }
                unit_count++;
            }
        }
		return unit_count;
	}else{
		while(!found){
			if(units[unit_count].type == type && our_unit_num == work_count){
				found=true;
			}else{
				if(units[unit_count].type != 'C'){
					 work_count++;
				}
				unit_count++;
			}
		}
		return unit_count;
	}
}

void send_num_data(struct unit units[], struct locators relative_positions[],
				   long long grid_size, double **p_variables_data,
				   double **p_variables_recv){
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int unit_count = find_unit(units, relative_positions, world_rank, 'B');
	int total_coupler_unit_count = units[unit_count].coupler_ranks.size();
    int ranks_per_coupler;
	double nodes_size;
	int buff_size = std::round(0.05 * (double) grid_size/8 * NVAR);
	*p_variables_data = (double*) malloc(buff_size * sizeof(double));
	*p_variables_recv = (double*) malloc(buff_size * sizeof(double));
	for(int i = 0; i < buff_size; i++){
		(*p_variables_data)[i] = std::rand()/1000;
		(*p_variables_recv)[i] = std::rand()/1000;
	}
    for(int z = 0; z < total_coupler_unit_count; z++){
		ranks_per_coupler = units[unit_count].coupler_ranks[z].size();
		int coupler_rank = units[unit_count].coupler_ranks[z][0];
		int unit_count_2 = find_unit(units, relative_positions, coupler_rank, 'C');
		if(units[unit_count_2].coupling_type == 'C'){
			nodes_size = std::round(0.018 * (double) grid_size/8);
		}else if(units[unit_count_2].coupling_type == 'O'){
			nodes_size = std::round(0.05 * (double) grid_size/8);
		}
		for(int z2 = 0; z2 < ranks_per_coupler; z2++){
			//Send the node sizes to each of the coupler ranks of each coupler unit.
			MPI_Send(&nodes_size, 1, MPI_DOUBLE, units[unit_count].coupler_ranks[z][z2], 0, MPI_COMM_WORLD);
		}
    }
}

void send_recv_data(struct unit units[], struct locators relative_positions[], 
					long long grid_size, int cycle_num, int total_cycles, 
					double *p_variables_data, double *p_variables_recv, FILE *fp){
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int unit_count = find_unit(units, relative_positions, world_rank, 'B');
	int total_coupler_unit_count = units[unit_count].coupler_ranks.size();
    fprintf(fp, "PIC cycle %d of %d comms starting\n", (cycle_num/pic_conversion_factor)+1, 
				total_cycles);
    for(int j = 0; j < total_coupler_unit_count; j++){
		//find the coupler unit index
		int coupler_rank = units[unit_count].coupler_ranks[j][0];
        int unit_count_2 = find_unit(units, relative_positions, coupler_rank, 'C');
		double send_size;
		//work out the amount of data to send
		if(units[unit_count_2].coupling_type == 'C'){
            send_size = std::round(0.018 * (double) grid_size/8);
        }else if(units[unit_count_2].coupling_type == 'O'){
			send_size = std::round(0.05 * (double) grid_size/8);
        }
		//send all the data to the couplers
        if(hide_search == true && pic_conversion_factor != 1){
			if((cycle_num % (search_freq*pic_conversion_factor)) == 0){
				MPI_Send(p_variables_data, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
			}else if((cycle_num % pic_conversion_factor) == pic_conversion_factor - 1){
				MPI_Recv(p_variables_recv, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}else{
				MPI_Send(p_variables_data, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
				MPI_Recv(p_variables_recv, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}else{
			MPI_Send(p_variables_data, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
			MPI_Recv(p_variables_recv, send_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
	fprintf(fp, "PIC cycle %d of %d comms ending\n", (cycle_num/pic_conversion_factor)+1, 
				total_cycles);
}

