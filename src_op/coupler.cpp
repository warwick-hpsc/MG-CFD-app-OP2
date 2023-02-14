#include "mgcfd_lib.h"
#ifdef deffenics 
    #include "dolfinx_lib.h"
#endif
#ifdef defsimpic
    #include "simpic_lib.h"
#endif
#include <stdio.h>
#include <mpi.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "../src/structures.h"
#include "../src/config.h"
#include "../src/const.h"
#include <time.h>
#include <iterator>
#include <iostream>
#include <chrono>
#include <map>

int main(int argc, char** argv){
	#include "coupler_config.h"
	FILE *ifp = fopen("cpx_input.cfg", "r");

	if(ifp == NULL){
		fprintf(stderr, "Can't open input file cpx_input.cfg\n");
		exit(1);
	}

	char coupler[] = "COUPLER";
	char mgcfd[] = "MG-CFD";
	char fenics[] = "FENICS";
	char simpic[] = "SIMPIC";
	char unit_1[] = "UNIT_1";
	char unit_2[] = "UNIT_2";
	char total[] = "TOTAL";
	char type[] = "TYPE";
	int mpi_ranks = 0;
	char keyword[8];//longest word is COUPLER
	char c_type[8]; //longest type is sliding
	int temp_unit;//temporarily stores how many processes each unit will have
	int num_of_units = 0;//number of coupler units as read from the TOTAL value
	int temp_count = 0;//used to count the total number of units and verify it matches the TOTAL value
	int coupler_count = 0;//used to count total number of coupler units
	int mgcfd_count = 0;//used to count total umber of MG-CFD units
	int fenics_count = 0;//used to count the number of FENICS units
	int simpic_count = 0;//used to count the number of CUP-CFD units

	fscanf(ifp, "%s %d", keyword, &temp_unit);
	if(strcmp(keyword, total) == 0){
		num_of_units = temp_unit;
	}

	struct unit units[num_of_units];//creates AoS for either coupler unit or MG-CFD unit

	while(fscanf(ifp, "%s %d", keyword, &temp_unit) != EOF){//filling the AoS with unit information
		if(strcmp(keyword, coupler) == 0){
			units[temp_count].type = 'C';
			units[temp_count].processes = temp_unit;
			fscanf(ifp, "%s %s", keyword, c_type);
			if(strcmp(keyword, type) != 0){
				fprintf(stderr, "Error: You must specify the type of a coupler after the definition, aborting... \n");
				exit(1);
			}
			if(strcmp(c_type, "SLIDING") == 0){
				units[temp_count].coupling_type = 'S';
			}else if(strcmp(c_type, "CHT") == 0){
				units[temp_count].coupling_type = 'C';
			}else if(strcmp(c_type, "OVERSET") == 0){
				units[temp_count].coupling_type = 'O';
			}else{
				fprintf(stderr, "Error: couplers must be of type CHT, OVERSET, or SLIDING, aborting... \n");
				exit(1);
			}
			temp_count++;
			coupler_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, mgcfd) == 0){
			units[temp_count].type = 'M';
			units[temp_count].processes = temp_unit;
			temp_count++;
			mgcfd_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, fenics) == 0){
			#ifndef deffenics
				fprintf(stderr, "Error: CPX has not been compiled with FEniCS X support.\n");
				exit(1);
			#endif				
			units[temp_count].type = 'F';
			units[temp_count].processes = temp_unit;
			temp_count++;
			fenics_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, simpic) == 0){
			#ifndef defsimpic
				fprintf(stderr, "Error: CPX has not been compiled with SIMPIC support.\n");
				exit(1);
			#endif
			units[temp_count].type = 'P';
			units[temp_count].processes = temp_unit;
			temp_count++;
			simpic_count++;
			mpi_ranks += temp_unit;
		}else if(strcmp(keyword, unit_1) == 0 || strcmp(keyword, unit_2) == 0){
			units[temp_count-1].mgcfd_units.push_back(temp_unit);
		}
	}
	if(temp_count != num_of_units){
		fprintf(stderr, "Error: there is a mismatch in the number of cpx/unit instances, aborting...\n");
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
		printf("\n Total number of units: %d\n No of MG-CFD units: %d\n No of FENICS units: %d\n No of SIMPIC units: %d\n No of Coupler units: %d\n\n Unit list:\n", num_of_units, mgcfd_count, fenics_count, simpic_count, coupler_count);

		int coupler_count_temp = 0;
		int work_count_temp = 0;

		for(int i = 0; i < num_of_units; i++){
			if(units[i].type == 'C'){
				coupler_count_temp++;
				printf("  Coupler unit %d\n  Coupler type %c\n  Number of ranks assigned: %d\n  with the two units %d and %d\n\n", coupler_count_temp, units[i].coupling_type, units[i].processes, units[i].mgcfd_units[0], units[i].mgcfd_units[1]);
			}else{
				work_count_temp++;
				printf("  Work unit %d type %c\n  Number of ranks assigned: %d\n\n", work_count_temp, units[i].type, units[i].processes);
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
	int temp_work = 0;

	std::vector<int> temp_coupler_position;
	std::vector<int> temp_work_position;

	for(int i=0; i<num_of_units;i++){
		if(units[i].type != 'C'){
			temp_work++;
			temp_work_position.push_back(temp_work);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = units[i].type;
				relative_positions[j].placelocator = temp_work;
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
				if(relative_positions[j].placelocator == (k2 + 1) && relative_positions[j].typelocator != 'C'){
                    temp_mgcfd_processes.push_back(j);//add mgcfd ranks to the temp_mgcfd_processes vector to transfer to the units[k & k2].coupler_ranks
				}
			}

			units[i].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
			units[k2].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 1st MG-CFD unit to the discovered ranks

			int k3 = units[i].mgcfd_units[1] - 1;//k3 is the 2nd MGCFD instance that the coupler unit manages

			units[k3].coupler_ranks.push_back(temp_coupler_processes);


			temp_mgcfd_processes.clear();
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == (k3 + 1) && relative_positions[j].typelocator != 'C'){
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
            if(units[i].coupler_ranks.size() != 0){
			    int no_of_elements_coupler = units[i].coupler_ranks[0].size();
			    for(int j = 0; j < no_of_elements_coupler; j++){
				    printf("Coupler ranks: %d\n", units[i].coupler_ranks[0][j]);
			    }
            }
			int no_of_elements = units[i].mgcfd_ranks[0].size();
			for(int j = 0; j < no_of_elements; j++){
			    printf("Unit Ranks: %d\n", units[i].mgcfd_ranks[0][j]);
			    if(units[i].type == 'C'){
				    printf("Unit Ranks: %d\n", units[i].mgcfd_ranks[1][j]);
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
	bool is_mgcfd = false;
	bool is_fenics = false;
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
		if(units[i].type != 'C'){
			int ranks[units[i].mgcfd_ranks[0].size()];
			std::copy(units[i].mgcfd_ranks[0].begin(), units[i].mgcfd_ranks[0].end(), ranks);
			MPI_Group_incl(world_group, units[i].mgcfd_ranks[0].size(), ranks, &new_groups[i]);
			MPI_Comm_create_group(MPI_COMM_WORLD, new_groups[i], 0, &new_comms[i]);
			for(int j=units[i].mgcfd_ranks[0][0]; j<units[i].mgcfd_ranks[0][0]+units[i].mgcfd_ranks[0].size();j++){//if rank of processes matches a rank for a particular MG-CFD unit, assign the new communicator
				if(rank == j){
					new_comm=new_comms[i];
					if(units[i].type == 'M'){
						is_mgcfd = true;
					} else if(units[i].type == 'F'){
						is_fenics = true;
					}
					instance_number = relative_positions[rank].placelocator;
				}
			}
		}
	}
    MPI_Fint comms_shell = MPI_Comm_c2f(new_comm);
	//end of the set up we then call mgcfd main or fenics main if its not a coupler.
	if(!is_coupler){
		if(is_mgcfd){
            main_mgcfd(argc, argv, comms_shell, instance_number, units, relative_positions);
		}else if(is_fenics){
			#ifdef deffenics
                main_dolfinx(argc, argv, comms_shell, instance_number, units, relative_positions);
			#endif
            MPI_Finalize();
		}else{
			#ifdef defsimpic
				//printf("launch simpic here");
				main_simpic(argc, argv, comms_shell, instance_number, units, relative_positions);
				//main_cup(argc, argv, comms_shell, instance_number, units, relative_positions);
            #endif
			MPI_Finalize();
		}
	}else{
		if(superdebug){
			#ifdef defsimpic
			#endif
			MPI_Finalize();
			
		}else{
			MPI_Comm coupler_comm = MPI_Comm_f2c(comms_shell);

			int my_rank;
	  		MPI_Comm_rank(coupler_comm, &my_rank);

			bool found = false;
			int unit_count = 0;
	    	while(!found){
				for(int j=units[unit_count].coupler_ranks[0][0]; j<units[unit_count].coupler_ranks[0][0]+units[unit_count].coupler_ranks[0].size();j++){//if rank of processes matches a rank for a particular coupler unit, assign the new communicator
					if(units[unit_count].type == 'C' && rank == j){
						found=true;
					}
				}
	            if(found == false){
	                unit_count++;
	            }
	    	}
			
			int left_rank = units[unit_count].mgcfd_ranks[0][0];
			int right_rank = units[unit_count].mgcfd_ranks[1][0];
			 
	        double left_nodes_size = 0.0;
	        double right_nodes_size = 0.0;
 
	        MPI_Recv(&left_nodes_size, 1, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(&right_nodes_size, 1, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int left_right_size = (int) ((left_nodes_size + right_nodes_size)/2);

			double adjusted_sizes_left = ceil((left_nodes_size/34000)/4);
			double adjusted_sizes_right = ceil((right_nodes_size/34000)/4);

			//set the scaling for the search and interpolation routines
			double left_search_scaling = 0;
			double right_search_scaling = 0;
			double search_repeats;
            double interp_scaling;
			if(ultrafastsearch){
				search_repeats = 0.25;
				interp_scaling = 90;
			}else{
				search_repeats = 10;
				interp_scaling = 3500;
			}
			
            int coupler_vars = 0;
            if(units[unit_count].coupling_type == 'S'){
                left_search_scaling = search_repeats * adjusted_sizes_left;
                right_search_scaling = search_repeats * adjusted_sizes_right;
                //interp_scaling = 3500/((adjusted_sizes_left + adjusted_sizes_right)/2);
                coupler_vars = 5;
            }else if(units[unit_count].coupling_type == 'O'){
				left_search_scaling = search_repeats * adjusted_sizes_left;
                right_search_scaling = search_repeats * adjusted_sizes_right;
                //interp_scaling = 3500/((adjusted_sizes_left + adjusted_sizes_right)/2);
                coupler_vars = 1;
			}else if(units[unit_count].coupling_type == 'C'){
                left_search_scaling = 3;
                right_search_scaling = 2;
                coupler_vars = 1;
            }

	        double *left_p_variables_recv = (double *) malloc(left_nodes_size * coupler_vars * sizeof(double));
	        double *right_p_variables_recv = (double *) malloc(right_nodes_size * coupler_vars * sizeof(double));
			double *left_p_variables = (double *) malloc(left_right_size * coupler_vars * sizeof(double));
			double *right_p_variables = (double *) malloc(left_right_size * coupler_vars * sizeof(double));

			int total_ranks = units[unit_count].coupler_ranks[0].size();
			int root_rank = units[unit_count].coupler_ranks[0][0];

			double left_right_size_chunks = left_right_size / total_ranks;

			double vector_counter;
			int sub_count;
			double vector_counter_max;

			//p_variables storage for scatter/gather
			double *left_p_variables_sg = (double *) malloc((left_right_size_chunks) * coupler_vars * sizeof(double));
	        double *right_p_variables_sg = (double *) malloc((left_right_size_chunks) * coupler_vars * sizeof(double));

			//p_variable vectors to represent brute force search
			std::vector< std::vector<double> > left_vector_of_state_vars;
			std::vector< std::vector<double> > right_vector_of_state_vars;

			//p_variable vectors to represent tree-based search
			std::map< int, std::vector<double> > left_map_of_state_vars;
			std::map< int, std::vector<double> > right_map_of_state_vars;

			std::vector<double> node_state_vars;
			std::vector<double> node_state_vars_left;
			std::vector<double> node_state_vars_right;
			std::vector<double> node_state_vars_temp;

			//set up some random data for cht interpolation
			int ar_size_max = left_right_size*0.9;
			if(rank == root_rank && debug == true){
				printf("size is %d\n", ar_size_max);
			}
			srand((unsigned)time(NULL));
			double data_ran[ar_size_max][4];
			for(int i = 0; i < ar_size_max; i++){
				for(int k = 0; k < 4; k++){
					data_ran[i][k] = rand()/1000;
				}
			}
			
			std::chrono::duration<double> total_seconds;
			std::chrono::duration<double> non_coupling_secs;
			std::chrono::duration<double> pure_compute_sec;
			std::chrono::duration<double> wait_sec;
			std::chrono::time_point<std::chrono::steady_clock> start;
			std::chrono::time_point<std::chrono::steady_clock> start1;

			for(int cycle_counter = 0; cycle_counter < coupler_cycles; cycle_counter++){
				int local_size;
				MPI_Comm_size(coupler_comm, &local_size);
				if(rank == root_rank){
					printf("Coupler cycle %d starting\n", cycle_counter+1);
					start = std::chrono::steady_clock::now();
					MPI_Recv(left_p_variables_recv, left_nodes_size * coupler_vars, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					start1 = std::chrono::steady_clock::now();
					MPI_Recv(right_p_variables_recv, right_nodes_size * coupler_vars, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					auto end = std::chrono::steady_clock::now();
					non_coupling_secs += (end-start);
					wait_sec = (end-start1);
					//sort the left and right variables to be roughly the avg of the two sides
					int counter = 0;
					for(int i = 0; i < left_nodes_size * coupler_vars; i++){
						if(counter < left_right_size * coupler_vars){
							left_p_variables[counter] = left_p_variables_recv[i];
						}else if(counter < 2 * left_nodes_size * coupler_vars){
							right_p_variables[counter-(left_right_size * coupler_vars)] = left_p_variables_recv[i];
						}
						counter++;
					}
					for(int i = 0; i < right_nodes_size * coupler_vars; i++){
						if(counter < left_right_size * coupler_vars){
							left_p_variables[counter] = right_p_variables_recv[i];
						}else if(counter < 2 * left_right_size * coupler_vars){
							right_p_variables[counter-(left_right_size * coupler_vars)] = right_p_variables_recv[i];
						}
						counter++;
					}
					start = std::chrono::steady_clock::now();
		        }

				MPI_Barrier(coupler_comm);
				MPI_Scatter(left_p_variables, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, left_p_variables_sg, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, 0, coupler_comm);
				MPI_Barrier(coupler_comm);
				MPI_Scatter(right_p_variables, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, right_p_variables_sg, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, 0, coupler_comm);
				
				if(MUM == 0){
					MPI_Bcast(left_p_variables, left_right_size * coupler_vars, MPI_DOUBLE, 0, coupler_comm);
					MPI_Bcast(right_p_variables, left_right_size * coupler_vars, MPI_DOUBLE, 0, coupler_comm);
				}

				if(rank == root_rank){
					start1 = std::chrono::steady_clock::now();
				}
				
				//rendezvous routines start
				if(units[unit_count].coupling_type == 'S' || cycle_counter == 0){
					if((cycle_counter % search_freq) == 0){
						if(MUM == 0){
							vector_counter_max = left_right_size;//this is size of mesh recieved from broadcast
						}else{
							vector_counter_max = left_right_size_chunks;//this is size of mesh recieved from scatter
						}
						if(!fastsearch){
							for(int l = 0; l < left_search_scaling; l++){
								vector_counter = 0;
								sub_count = 0;
								while(sub_count < total_ranks){
									left_vector_of_state_vars.clear();
									while(vector_counter < (vector_counter_max/total_ranks)){
										node_state_vars.clear();
										for(int i = 0; i<coupler_vars; i++){
											if(MUM == 0){
												node_state_vars.push_back(*(left_p_variables + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along left_p_variables in chunks of NVAR
											}else{
												node_state_vars.push_back(*(left_p_variables_sg + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along left_p_variables in chunks of NVAR
											}
										}
										left_vector_of_state_vars.insert(left_vector_of_state_vars.begin(), node_state_vars);
										vector_counter++;
									}
									vector_counter = 0;
									sub_count++;
								}
							}
						}else{
							for(int l = 0; l < left_search_scaling; l++){
								left_map_of_state_vars.clear();
								vector_counter = 0;
								while(vector_counter < (vector_counter_max)){
									node_state_vars.clear();
									for(int i = 0; i<coupler_vars; i++){
										if(MUM == 0){
											node_state_vars.push_back(*(left_p_variables + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along left_p_variables in chunks of NVAR
										}else{
											node_state_vars.push_back(*(left_p_variables_sg + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along left_p_variables in chunks of NVAR
										}
									}
									left_map_of_state_vars.insert(std::make_pair(vector_counter, node_state_vars));
									vector_counter++;
								}
							}
				        }
					}

					if((cycle_counter % search_freq) == 0){
						if(MUM == 0){
							vector_counter_max = left_right_size;//this is size of mesh recieved from broadcast
						}else{
							vector_counter_max = left_right_size_chunks;//this is size of mesh recieved from scatter
						}
						if(!fastsearch){
							for(int l = 0; l < right_search_scaling; l++){
								vector_counter = 0;
								sub_count = 0;
								while(sub_count < total_ranks){
									right_vector_of_state_vars.clear();
									while(vector_counter < (vector_counter_max/total_ranks)){
										node_state_vars.clear();
										for(int i = 0; i<coupler_vars; i++){
											if(MUM == 0){
												node_state_vars.push_back(*(right_p_variables + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along right_p_variables in chunks of NVAR
											}else{
												node_state_vars.push_back(*(right_p_variables_sg + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along right_p_variables in chunks of NVAR
											}
										}
										right_vector_of_state_vars.insert(right_vector_of_state_vars.begin(), node_state_vars);
										vector_counter++;
									}
									vector_counter = 0;
									sub_count++;
								}
							}
						}else{
							for(int l = 0; l < right_search_scaling; l++){
								right_map_of_state_vars.clear();
								vector_counter = 0;
								while(vector_counter < (vector_counter_max)){
									node_state_vars.clear();
									for(int i = 0; i<coupler_vars; i++){
										if(MUM == 0){
											node_state_vars.push_back(*(right_p_variables + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along right_p_variables in chunks of NVAR
										}else{
											node_state_vars.push_back(*(right_p_variables_sg + (static_cast<long long>(vector_counter) * coupler_vars) + i));//essentially move along right_p_variables in chunks of NVAR
										}
									}
									right_map_of_state_vars.insert(std::make_pair(vector_counter, node_state_vars));
									vector_counter++;
								}
							}
				        }
					}
					//this routine is used to convert the tree in the based search to a vector so the interpolation routine can run as before
					if(fastsearch){
						vector_counter_max = left_right_size_chunks;
						left_vector_of_state_vars.clear();
						right_vector_of_state_vars.clear();
						for(int i=0; i<(vector_counter_max/total_ranks); i++){
							left_vector_of_state_vars.insert(left_vector_of_state_vars.end(), left_map_of_state_vars.at(i));
							right_vector_of_state_vars.insert(right_vector_of_state_vars.end(), right_map_of_state_vars.at(i));
						}
						if(units[unit_count].coupling_type == 'S'){
							left_search_scaling = adjusted_sizes_left;
							right_search_scaling = adjusted_sizes_right;
						}
					}
				}
				//rendezvous routines end
	
				//interpolate routine start
				if(units[unit_count].coupling_type == 'S' || units[unit_count].coupling_type == 'O'){
					node_state_vars_left.clear();
					node_state_vars_right.clear();
					node_state_vars_temp.clear();
					vector_counter = 0;
					sub_count = 0;
					vector_counter_max = left_right_size_chunks;//this is size of mesh recieved from scatter
					while(sub_count < (total_ranks*(interp_scaling/((adjusted_sizes_left + adjusted_sizes_right)/2)))){//TODO: changes the adjusted_sizes to whichever is lower
						while(vector_counter < vector_counter_max/total_ranks){
							node_state_vars_left = left_vector_of_state_vars.at(vector_counter);
							node_state_vars_right = right_vector_of_state_vars.at(vector_counter);
							node_state_vars_temp = node_state_vars_right;
							for(int i = 0; i<coupler_vars; i++){
								node_state_vars_right.at(i) = (node_state_vars_left.at(i) + node_state_vars_right.at(i))/2;
								node_state_vars_left.at(i) = (node_state_vars_left.at(i) + node_state_vars_temp.at(i))/2;
							}
							vector_counter++;
						}
						vector_counter = 0;
						sub_count++;
					}
				}else if(units[unit_count].coupling_type == 'C'){
					std::vector<double> quad_array_rt;
					std::vector<double> quad_array_lt;
					for(int l = 0; l < 1; l++){
						quad_array_rt.clear();
						quad_array_lt.clear();
						for(int i = 0; i < ar_size_max*0.7; i++){
							quad_array_rt.push_back(0.0);
							quad_array_lt.push_back(0.0);
							for(int k = 0; k < 4; k++){
								quad_array_rt.at(i) = (quad_array_rt.at(i) + data_ran[i][k]);
								quad_array_lt.at(i) = (quad_array_lt.at(i) + data_ran[i][k]);
							}
							quad_array_rt.at(i) = (quad_array_rt.at(i)/4);
							quad_array_lt.at(i) = (quad_array_lt.at(i)/4);
						}
					}
					vector_counter_max = left_right_size_chunks;//this is size of mesh recieved from scatter
					for(int l = 0; l < 7; l++){
						for(int b = 0; b < vector_counter_max/5; b++){
							int temp = 0;
							for(int i = 0; i < 3; i++){
								if(data_ran[b][i] != data_ran[b][i % 3 + 1]){
									temp = temp + 1;
									quad_array_rt.at(b % (std::max(ar_size_max/20,1))) = data_ran[b][temp];
								}
							}
						}
					}
				}
				if(rank == root_rank){
					auto end1 = std::chrono::steady_clock::now();
					pure_compute_sec = (end1-start1);
				}
				
				//interpolate routine end
				MPI_Barrier(coupler_comm);
		        MPI_Gather(left_p_variables_sg, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, left_p_variables, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, 0, coupler_comm);
				MPI_Barrier(coupler_comm);
		        MPI_Gather(right_p_variables_sg, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, right_p_variables, (left_right_size_chunks * coupler_vars), MPI_DOUBLE, 0, coupler_comm);
				
				if(rank == root_rank){
		            MPI_Send(right_p_variables_recv, right_nodes_size * coupler_vars, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
		            MPI_Send(left_p_variables_recv, left_nodes_size * coupler_vars, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
					auto end = std::chrono::steady_clock::now();
					total_seconds += (end-start);
					printf("Coupler cycle %d ending\n", cycle_counter+1);
		        }
			}
			MPI_Barrier(coupler_comm);
			if(rank == root_rank){
				printf("total coupling time is %f\n", total_seconds.count());
				printf("total time waiting %f\n", non_coupling_secs.count());
				printf("time between left send and right send %f\n",wait_sec.count());
				printf("total pure compute time is %f\n", pure_compute_sec.count());
			}
			MPI_Finalize();
	   		exit(0);
		}
	}
}

