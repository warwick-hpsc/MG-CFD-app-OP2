#include "mgcfd_lib.h"
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "../src/structures.h"
#include "../src/config.h"


int main(int argc, char** argv){

	FILE *ifp = fopen("cpx_input.cfg", "r");

	if(ifp == NULL){
		fprintf(stderr, "Can't open input file cpx_input.cfg\n");
		exit(1);
	}

	char coupler[] = "COUPLER";
	char mgcfd[] = "MG-CFD";
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
		if(units[i].type == 'C'){
			temp_coupler++;
			temp_coupler_position.push_back(temp_coupler);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = 'C';
				relative_positions[j].placelocator = temp_coupler;
			}
			temp_marker += units[i].processes;
		}
		if(units[i].type == 'M'){
			temp_mgcfd++;
			temp_mgcfd_position.push_back(temp_mgcfd);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = 'M';
				relative_positions[j].placelocator = temp_mgcfd;
			}
			temp_marker += units[i].processes;	

			if(temp_coupler_position.size() > 0 && temp_mgcfd_position.size() == 2){//if there is an unassigned coupler unit and 2 unassigned mg-cfd units
				int k = 0;
				while(units[k].coupler_ranks.empty() == false || units[k].type == 'M'){//find the first unassigned coupler unit
					k++;
				}

				std::vector<int> temp_coupler_processes;//stores the process list to go into the units[k].coupler_ranks vector
				int temp_coupler_relative_position = temp_coupler_position.back();//the placelocator for the coupler unit
				temp_coupler_position.pop_back();
				for(int j=0; j<temp_marker; j++){
					if(relative_positions[j].placelocator == temp_coupler_relative_position && relative_positions[j].typelocator == 'C'){
						temp_coupler_processes.push_back(j);//add coupler ranks to the temp_coupler_processes vector to transfer to the units[k].coupler_ranks
					}
				}
				units[k].coupler_ranks = temp_coupler_processes;

				int k2 = 0;//keep k to store the mgcfd processes in the coupler unit 'mgcfd_ranks' vector
				while(units[k2].coupler_ranks.empty() == false || units[k2].type == 'C'){//find the first unassigned mgcfd unit
					k2++;
				}

				units[k2].coupler_ranks = temp_coupler_processes;

				std::vector<int> temp_mgcfd_processes;
				int temp_mgcfd_relative_position = temp_mgcfd_position.back();//the placelocator for the 1st MG-CFD unit
				temp_mgcfd_position.pop_back();
				for(int j=0; j<temp_marker; j++){
					if(relative_positions[j].placelocator == temp_mgcfd_relative_position && relative_positions[j].typelocator == 'M'){
						temp_mgcfd_processes.push_back(j);//add mgcfd ranks to the temp_mgcfd_processes vector to transfer to the units[k & k2].coupler_ranks
					}
				}

				units[k].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
				units[k2].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 1st MG-CFD unit to the discovered ranks

				int k3 = 0;//k3 is the 2nd MGCFD instance that the coupler unit manages
				while(units[k3].coupler_ranks.empty() == false || units[k3].type == 'C'){//find the second unassigned mgcfd unit
					k3++;
				}

				units[k3].coupler_ranks = temp_coupler_processes;
				temp_mgcfd_processes.clear();
				temp_mgcfd_relative_position = temp_mgcfd_position.back();//the placelocator for the 2nd MG-CFD unit
				temp_mgcfd_position.pop_back();
				for(int j=0; j<temp_marker; j++){
					if(relative_positions[j].placelocator == temp_mgcfd_relative_position && relative_positions[j].typelocator == 'M'){
						temp_mgcfd_processes.push_back(j);//add coupler ranks to the temp_coupler_processes vector to transfer to the units[k & k3].coupler_ranks
					}
				}

				units[k].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
				units[k3].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 2nd MG-CFD unit to the discovered ranks
			}
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
			printf("Coupler rank: %d\n", units[i].coupler_ranks[0]);
			int no_of_elements = units[i].mgcfd_ranks[0].size();
			for(int j = 0; j < no_of_elements; j++){
				printf("MGCFD Ranks: %d\n", units[i].mgcfd_ranks[0][j]);
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
			int ranks[1] = {units[i].coupler_ranks[0]};
			MPI_Group_incl(world_group, 1, ranks, &new_groups[i]);
			MPI_Comm_create_group(MPI_COMM_WORLD, new_groups[i], 0, &new_comms[i]);
			if(rank == ranks[0]){
				new_comm=new_comms[i];
				is_coupler = true;
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
		MPI_Comm coupler_comm = MPI_Comm_f2c(comms_shell);

		int cycle_counter = 0;

		bool found = false;
		int unit_count = 0;
    	while(!found){ /* Currently this is just for one coupler unit, so stop when the unit type is C*/
			if(units[unit_count].type == 'C'){
				found=true;
			}else{
				unit_count++;
			}
    	}
		int left_rank = units[unit_count].mgcfd_ranks[0][0];
		int left_size = static_cast<int>(units[unit_count].mgcfd_ranks[0].size());
		int right_rank = units[unit_count].mgcfd_ranks[1][0];
		int right_size = static_cast<int>(units[unit_count].mgcfd_ranks[1].size());
		int *left_rank_storage = new int[left_size];
		int *right_rank_storage = new int[right_size];

		while(cycle_counter < 25){/* TODO: get initial MPI message to get number of cycles */
			MPI_Recv(left_rank_storage, left_size, MPI_INT, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);/* currently ignore status - can be changed */
			/* Interpolate logic goes here */
			MPI_Send(left_rank_storage, left_size, MPI_INT, right_rank, 0, MPI_COMM_WORLD);

			MPI_Recv(right_rank_storage, right_size, MPI_INT, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);/* currently ignore status - can be changed */
			/* Interpolate logic goes here */
			MPI_Send(right_rank_storage, right_size, MPI_INT, left_rank, 0, MPI_COMM_WORLD);

			cycle_counter++;

		}
		MPI_Finalize();
   		exit(0);
	}
    
}

