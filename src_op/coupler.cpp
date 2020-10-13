#include "mgcfd_lib.h"
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "../src/structures.h"
#include "../src/config.h"
#include "../src/const.h"



int main(int argc, char** argv){

	FILE *ifp = fopen("cpx_input.cfg", "r");

	if(ifp == NULL){
		fprintf(stderr, "Can't open input file cpx_input.cfg\n");
		exit(1);
	}

	char coupler[] = "COUPLER";
	char mgcfd[] = "MG-CFD";
	char unit_1[] = "UNIT_1";
	char unit_2[] = "UNIT_2";
	char total[] = "TOTAL";
	bool debug = false;

	int mpi_ranks = 0;
	char keyword[8];//longest word is COUPLER
	int temp_unit;//temporarily stores how many processes each unit will have
	int num_of_units = 0;//number of coupler units as read from the TOTAL value
	int temp_count = 0;//used to count the total number of units and verify it matches the TOTAL value
	int coupler_count = 0;//used to count total number of coupler units
	int mgcfd_count = 0;//used to count total umber of MG-CFD units

	fscanf(ifp, "%s %d", keyword, &temp_unit);
	if(strcmp(keyword, total) == 0){
		num_of_units = temp_unit;
	}

	struct unit units[num_of_units];//creates AoS for either coupler unit or MG-CFD unit

	while(fscanf(ifp, "%s %d", keyword, &temp_unit) != EOF){//filling the AoS with unit information
		if(strcmp(keyword, coupler) == 0){
			units[temp_count].type = 'C';
			units[temp_count].processes = temp_unit;
			temp_count++;
			coupler_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, mgcfd) == 0){
			units[temp_count].type = 'M';
			units[temp_count].processes = temp_unit;
			temp_count++;
			mgcfd_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, unit_1) == 0 || strcmp(keyword, unit_2) == 0){
			units[temp_count-1].mgcfd_units.push_back(temp_unit);
		}
	}

	if(temp_count != num_of_units){
		fprintf(stderr, "Error: there is a mismatch in the number of cpx/mg-cfd instances, aborting...\n");
		exit(1);
	}

	fclose(ifp);
    MPI_Init(&argc, &argv);

    int rank, size;

    //get initial ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    if(rank == 0){
        printf("It's coupler time ;)");

		printf("\n Total number of units: %d\n No of MG-CFD units: %d\n No of Coupler units: %d\n\n Unit list:\n", num_of_units, mgcfd_count,coupler_count);

		int coupler_count_temp = 0;
		int mgcfd_count_temp = 0;

		for(int i = 0; i < num_of_units; i++){
			if(units[i].type == 'C'){
				coupler_count_temp++;
				printf("  Coupler unit %d\n  Number of ranks assigned: %d\n\n", coupler_count_temp, units[i].processes);
			}else{
				mgcfd_count_temp++;
				printf("  MG-CFD unit %d\n  Number of ranks assigned: %d\n\n", mgcfd_count_temp, units[i].processes);
			}
		}

        if(mpi_ranks != size){
            fprintf(stderr, "Error: there is a mismatch in the number requested MPI ranks and the number of ranks requested in the config file.\n");
            exit(1);
	    }
    }

	struct locators relative_positions[mpi_ranks];//creates AoS with an entry for every process

	int temp_coupler = 0;
	int temp_marker = 0;//to mark where in the relative_positions array we are
	int temp_mgcfd = 0;

	std::vector<int> temp_coupler_position;
	std::vector<int> temp_mgcfd_position;

	for(int i=0; i<num_of_units;i++){
		if(units[i].type == 'M'){
			temp_mgcfd++;
			temp_mgcfd_position.push_back(temp_mgcfd);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = 'M';
				relative_positions[j].placelocator = temp_mgcfd;
			}
			temp_marker += units[i].processes;	
		}
		if(units[i].type == 'C'){
			temp_coupler++;
			temp_coupler_position.push_back(temp_coupler);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = 'C';
				relative_positions[j].placelocator = temp_coupler;
			}
			temp_marker += units[i].processes;	

			std::vector<int> temp_coupler_processes;//stores the process list to go into the units[k].coupler_ranks vector
			int temp_coupler_relative_position = temp_coupler_position.back();//the placelocator for the coupler unit
			temp_coupler_position.pop_back();
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == temp_coupler_relative_position && relative_positions[j].typelocator == 'C'){
					temp_coupler_processes.push_back(j);//add coupler ranks to the temp_coupler_processes vector to transfer to the units[k].coupler_ranks
				}
			}
			units[i].coupler_ranks.push_back(temp_coupler_processes);

			int k2 = units[i].mgcfd_units[0] - 1;//keep k to store the mgcfd processes in the coupler unit 'mgcfd_ranks' vector

        
			units[k2].coupler_ranks.push_back(temp_coupler_processes);
            

			std::vector<int> temp_mgcfd_processes;
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == (k2 + 1) && relative_positions[j].typelocator == 'M'){
					temp_mgcfd_processes.push_back(j);//add mgcfd ranks to the temp_mgcfd_processes vector to transfer to the units[k & k2].coupler_ranks
				}
			}

			units[i].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
			units[k2].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 1st MG-CFD unit to the discovered ranks

			int k3 = units[i].mgcfd_units[1] - 1;//k3 is the 2nd MGCFD instance that the coupler unit manages


			units[k3].coupler_ranks.push_back(temp_coupler_processes);


			temp_mgcfd_processes.clear();
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == (k3 + 1) && relative_positions[j].typelocator == 'M'){
					temp_mgcfd_processes.push_back(j);//add coupler ranks to the temp_coupler_processes vector to transfer to the units[k & k3].coupler_ranks
				}
			}

			units[i].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
			units[k3].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 2nd MG-CFD unit to the discovered ranks
		} 
	}
	//for debugging purposes
	if(rank == 0 && debug == true){
		for(int i = 0; i<mpi_ranks; i++){
			printf("\n\nRank %d\n", i);
			printf("Type: %c\n",relative_positions[i].typelocator);
			printf("Place: %d\n",relative_positions[i].placelocator);
		}

		for(int i = 0; i<num_of_units; i++){
			printf("\n\nType %c\n", units[i].type);
			int no_of_elements_coupler = units[i].coupler_ranks[0].size();
			for(int j = 0; j < no_of_elements_coupler; j++){
				printf("Coupler ranks: %d\n", units[i].coupler_ranks[0][j]);
			}
			int no_of_elements = units[i].mgcfd_ranks[0].size();
			for(int j = 0; j < no_of_elements; j++){
				printf("MGCFD Ranks: %d\n", units[i].mgcfd_ranks[0][j]);
				if(units[i].type == 'C'){
					printf("MGCFD Ranks: %d\n", units[i].mgcfd_ranks[1][j]);
				}
			}
		}
	}
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	MPI_Group new_groups[num_of_units];
	MPI_Comm new_comms[num_of_units];
	MPI_Comm new_comm;
	bool is_coupler = false;
	int instance_number = 0;

	for(int i=0; i<num_of_units;i++){
		if(units[i].type == 'C'){
			int ranks[units[i].coupler_ranks[0].size()];
			std::copy(units[i].coupler_ranks[0].begin(), units[i].coupler_ranks[0].end(), ranks);
			MPI_Group_incl(world_group, units[i].coupler_ranks[0].size(), ranks, &new_groups[i]);
			MPI_Comm_create_group(MPI_COMM_WORLD, new_groups[i], 0, &new_comms[i]);
			for(int j=units[i].coupler_ranks[0][0]; j<units[i].coupler_ranks[0][0]+units[i].coupler_ranks[0].size();j++){//if rank of processes matches a rank for a particular coupler unit, assign the new communicator
				if(rank == j){
					new_comm=new_comms[i];
					is_coupler = true;
					instance_number = relative_positions[rank].placelocator;
				}
			}
		}
		if(units[i].type == 'M'){
			int ranks[units[i].mgcfd_ranks[0].size()];
			std::copy(units[i].mgcfd_ranks[0].begin(), units[i].mgcfd_ranks[0].end(), ranks);
			MPI_Group_incl(world_group, units[i].mgcfd_ranks[0].size(), ranks, &new_groups[i]);
			MPI_Comm_create_group(MPI_COMM_WORLD, new_groups[i], 0, &new_comms[i]);
			for(int j=units[i].mgcfd_ranks[0][0]; j<units[i].mgcfd_ranks[0][0]+units[i].mgcfd_ranks[0].size();j++){//if rank of processes matches a rank for a particular MG-CFD unit, assign the new communicator
				if(rank == j){
					new_comm=new_comms[i];
					instance_number = relative_positions[rank].placelocator;
				}
			}
		}
	}
    MPI_Fint comms_shell = MPI_Comm_c2f(new_comm);
	if(!is_coupler){
		main_mgcfd(argc, argv, comms_shell, instance_number, units, relative_positions);
	}else{
		#include "coupler_config.h"
		MPI_Comm coupler_comm = MPI_Comm_f2c(comms_shell);
		int cycle_counter = 0;

		bool found = false;
		int unit_count = 0;
    	while(!found){
			for(int j=units[unit_count].coupler_ranks[0][0]; j<units[unit_count].coupler_ranks[0][0]+units[unit_count].coupler_ranks[0].size();j++){//if rank of processes matches a rank for a particular coupler unit, assign the new communicator
				if(units[unit_count].type == 'C' && rank == j){
					found=true;
				}else if(!(j < (units[unit_count].coupler_ranks[0][0]+units[unit_count].coupler_ranks[0].size() - 1)) && found == false){//to make sure all ranks of each coupler units are found
					unit_count++;
				}
			}
    	}

		int left_rank = units[unit_count].mgcfd_ranks[0][0];
		int left_size = static_cast<int>(units[unit_count].mgcfd_ranks[0].size());
		int right_rank = units[unit_count].mgcfd_ranks[1][0];
		int right_size = static_cast<int>(units[unit_count].mgcfd_ranks[1].size());
		int *left_rank_storage = new int[left_size];
		int *right_rank_storage = new int[right_size];
        
        double left_nodes_sizes[4];
        double right_nodes_sizes[4];
        MPI_Recv(left_nodes_sizes, 4, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(right_nodes_sizes, 4, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double *left_p_variables_l0, *left_p_variables_l1, *left_p_variables_l2, *left_p_variables_l3;
        double *right_p_variables_l0, *right_p_variables_l1, *right_p_variables_l2, *right_p_variables_l3;

        left_p_variables_l0 = (double *) malloc(left_nodes_sizes[0] * NVAR * sizeof(double));
        left_p_variables_l1 = (double *) malloc(left_nodes_sizes[1] * NVAR * sizeof(double));
        left_p_variables_l2 = (double *) malloc(left_nodes_sizes[2] * NVAR * sizeof(double));
        left_p_variables_l3 = (double *) malloc(left_nodes_sizes[3] * NVAR * sizeof(double));

        right_p_variables_l0 = (double *) malloc(right_nodes_sizes[0] * NVAR * sizeof(double));
        right_p_variables_l1 = (double *) malloc(right_nodes_sizes[1] * NVAR * sizeof(double));
        right_p_variables_l2 = (double *) malloc(right_nodes_sizes[2] * NVAR * sizeof(double));
        right_p_variables_l3 = (double *) malloc(right_nodes_sizes[3] * NVAR * sizeof(double));

		/* Calculate sizes here */


		int total_ranks = units[unit_count].coupler_ranks[0].size();

    	double left_nodes_size_chunks[4];
		double left_extra_chunk[4];
		double right_nodes_size_chunks[4];
		double right_extra_chunk[4];
		double left_add_amount[4];
		double right_add_amount[4];
		bool chunk_add = true;

		for(int j = 0; j < 4; j++){
				left_nodes_size_chunks[j] = left_nodes_sizes[j] / total_ranks;
				left_extra_chunk[j] = left_nodes_sizes[j] - left_nodes_size_chunks[j]*total_ranks;//since dividing integers is a floor function in C++, there may be a little bit left which we add onto the first transfer
				right_nodes_size_chunks[j] = right_nodes_sizes[j] / total_ranks;
				right_extra_chunk[j] = right_nodes_sizes[j] - right_nodes_size_chunks[j]*total_ranks;
            }
		while((cycle_counter < 25) && ((cycle_counter % upd_freq) == 0)){/* Change this value to the number of cycles if it is not the default*/

			for(int k = 0; k < 4; k++){
				if(rank == units[unit_count].coupler_ranks[0][0]){//since the first rank of each coupler unit gets the extra chunk, check to see if the coupler unit process is the first rank
					left_add_amount[k] = left_nodes_size_chunks[k] + left_extra_chunk[k];
					right_add_amount[k] = right_nodes_size_chunks[k] + right_extra_chunk[k];
				}else{
					left_add_amount[k] = left_nodes_size_chunks[k];
					right_add_amount[k] = right_nodes_size_chunks[k];
				}
			}

            MPI_Recv(left_p_variables_l0, left_add_amount[0], MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(left_p_variables_l1, left_add_amount[1], MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(left_p_variables_l2, left_add_amount[2], MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(left_p_variables_l3, left_add_amount[3], MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(left_p_variables_l0, left_add_amount[0], MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
            MPI_Send(left_p_variables_l1, left_add_amount[1], MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
            MPI_Send(left_p_variables_l2, left_add_amount[2], MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
            MPI_Send(left_p_variables_l3, left_add_amount[3], MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);

            MPI_Recv(right_p_variables_l0, right_add_amount[0], MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(right_p_variables_l1, right_add_amount[1], MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(right_p_variables_l2, right_add_amount[2], MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(right_p_variables_l3, right_add_amount[3], MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(right_p_variables_l0, right_add_amount[0], MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
            MPI_Send(right_p_variables_l1, right_add_amount[1], MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
            MPI_Send(right_p_variables_l2, right_add_amount[2], MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
            MPI_Send(right_p_variables_l3, right_add_amount[3], MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);

			cycle_counter = cycle_counter + upd_freq;

		}
		MPI_Finalize();
   		exit(0);
	}
    
}




