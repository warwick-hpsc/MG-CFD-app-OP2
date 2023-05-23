#include "mgcfd_lib.h"
#ifdef deffenics 
    #include "dolfinx_lib.h"
#endif
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "../src/structures.h"
#include "../src/config.h"
#include "../src/const.h"
#include <time.h>
#include <iterator>
#include <iostream>
#include <map>

#define ASIO_STANDALONE
#include <asio.hpp>


int main(int argc, char** argv){

	FILE *ifp = fopen("cpx_input.cfg", "r");

	if(ifp == NULL){
		fprintf(stderr, "Can't open input file cpx_input.cfg\n");
		exit(1);
	}

	char coupler[] = "COUPLER";
	char mgcfd[] = "MG-CFD";
	char fenics[] = "FENICS";
	char unit_1[] = "UNIT_1";
	char unit_2[] = "UNIT_2";
	char total[] = "TOTAL";
	char type[] = "TYPE";
	char l_mgcfd[] = "L_MG-CFD";
	char l_coupler[] = "L_COUPLER";
	char l_unit[] = "L_UNIT";

	bool debug = true;

	int mpi_ranks = 0;
	char keyword[10];//longest word is COUPLER
	char c_type[8]; //longest type is sliding
	int temp_unit;//temporarily stores how many processes each unit will have
	int num_of_units = 0;//number of coupler units as read from the TOTAL value
	int temp_count = 0;//used to count the total number of units and verify it matches the TOTAL value
	int coupler_count = 0;//used to count total number of coupler units
	int mgcfd_count = 0;//used to count total umber of MG-CFD units
	int fenics_count = 0;//used to count the number of FENICS units
	int l_mgcfd_unit = -1;//used to find the loosely-coupled MG-CFD unit in the units array

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
			}else{
				fprintf(stderr, "Error: couplers must be of type CHT or SLIDING, aborting... \n");
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
		}else if(strcmp(keyword, unit_1) == 0 || strcmp(keyword, unit_2) == 0){
			units[temp_count-1].mgcfd_units.push_back(temp_unit);
		}else if(strcmp(keyword, l_coupler) == 0){
			units[temp_count].type = 'L';
			units[temp_count].processes = temp_unit;
			units[temp_count].coupling_type = 'S';
			temp_count++;
			coupler_count++;
			mpi_ranks += temp_unit;
		}else if (strcmp(keyword, l_unit) == 0){
			units[temp_count-1].mgcfd_units.push_back(temp_unit);
			//units[temp_count-1].mgcfd_units.push_back(-1);
		}else if(strcmp(keyword, l_mgcfd) == 0){
			l_mgcfd_unit = temp_unit;
		}
		
	}
	if(temp_count != num_of_units){
		fprintf(stderr, "Error: there is a mismatch in the number of cpx/mg-cfd instances, aborting...\n");
		fprintf(stderr, "temp_count = %d, num_of_units = %d\n", temp_count, num_of_units);
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
		printf("\n Total number of units: %d\n No of MG-CFD units: %d\n No of FENICS units: %d\n No of Coupler units: %d\n\n Unit list:\n", num_of_units, mgcfd_count, fenics_count, coupler_count);

		int coupler_count_temp = 0;
		int work_count_temp = 0;

		for(int i = 0; i < num_of_units; i++){
			if(units[i].type == 'C'){
				coupler_count_temp++;
				printf("  Coupler unit %d\n  Coupler type %c\n  Number of ranks assigned: %d\n  with the two units %d and %d\n\n", coupler_count_temp, units[i].coupling_type, units[i].processes, units[i].mgcfd_units[0], units[i].mgcfd_units[1]);
			}else if(units[i].type == 'L'){
				coupler_count_temp++;
				printf("  Loose Coupler unit %d\n  Number of ranks assigned: %d\n  with the unit %d\n\n", coupler_count_temp, units[i].processes, units[i].mgcfd_units[0]);
			}
			else{
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
		if(units[i].type == 'M' || units[i].type == 'F'){
			temp_work++;
			temp_work_position.push_back(temp_work);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = units[i].type;
				relative_positions[j].placelocator = temp_work;
			}
			temp_marker += units[i].processes;
		}       
		if(units[i].type == 'C' || units[i].type == 'L'){
			temp_coupler++;
			temp_coupler_position.push_back(temp_coupler);
			for(int j=temp_marker; j < temp_marker + units[i].processes; j++){
				relative_positions[j].typelocator = units[i].type;
				relative_positions[j].placelocator = temp_coupler;
			}
			temp_marker += units[i].processes;

			std::vector<int> temp_coupler_processes;//stores the process list to go into the units[k].coupler_ranks vector
			int temp_coupler_relative_position = temp_coupler_position.back();//the placelocator for the coupler unit
			temp_coupler_position.pop_back();
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == temp_coupler_relative_position && relative_positions[j].typelocator == units[i].type){
					temp_coupler_processes.push_back(j);//add coupler ranks to the temp_coupler_processes vector to transfer to the units[k].coupler_ranks
				}
			}
			units[i].coupler_ranks.push_back(temp_coupler_processes);

			int k2 = units[i].mgcfd_units[0] - 1;//keep k to store the mgcfd processes in the coupler unit 'mgcfd_ranks' vector


			units[k2].coupler_ranks.push_back(temp_coupler_processes);


			std::vector<int> temp_mgcfd_processes;
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == (k2 + 1) && relative_positions[j].typelocator != units[i].type){
                    temp_mgcfd_processes.push_back(j);//add mgcfd ranks to the temp_mgcfd_processes vector to transfer to the units[k & k2].coupler_ranks
				}
			}

			units[i].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the coupler unit to the discovered ranks
			units[k2].mgcfd_ranks.push_back(temp_mgcfd_processes);//set mgcfd_ranks of the 1st MG-CFD unit to the discovered ranks

			if(units[i].type == 'C'){
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

		if (l_mgcfd_unit > -1 && temp_work == l_mgcfd_unit && units[i].type == 'M') {
			std::vector<int> temp_mgcfd_processes;
			for(int j=0; j<temp_marker; j++){
				if(relative_positions[j].placelocator == l_mgcfd_unit){
                    temp_mgcfd_processes.push_back(j);//add mgcfd ranks to the temp_mgcfd_processes vector to transfer to the units[k & k2].coupler_ranks
				}
			}

			units[i].mgcfd_ranks.push_back(temp_mgcfd_processes);
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
	bool is_mgcfd = false;
	int instance_number = 0;

	for(int i=0; i<num_of_units;i++){
		if(units[i].type == 'C' || units[i].type == 'L'){
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
		if(units[i].type == 'M' || units[i].type == 'F'){
			int ranks[units[i].mgcfd_ranks[0].size()];
			std::copy(units[i].mgcfd_ranks[0].begin(), units[i].mgcfd_ranks[0].end(), ranks);
			MPI_Group_incl(world_group, units[i].mgcfd_ranks[0].size(), ranks, &new_groups[i]);
			MPI_Comm_create_group(MPI_COMM_WORLD, new_groups[i], 0, &new_comms[i]);
			for(int j=units[i].mgcfd_ranks[0][0]; j<units[i].mgcfd_ranks[0][0]+units[i].mgcfd_ranks[0].size();j++){//if rank of processes matches a rank for a particular MG-CFD unit, assign the new communicator
				if(rank == j){
					new_comm=new_comms[i];
					if(units[i].type == 'M'){
						is_mgcfd = true;
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
            main_mgcfd(argc, argv, comms_shell, instance_number, units, relative_positions, l_mgcfd_unit);
		}else{
            #ifdef deffenics
                main_dolfinx(argc, argv, comms_shell, instance_number, units, relative_positions);
			#endif
            MPI_Finalize();
		}
	}else{
		#include "coupler_config.h"
		MPI_Comm coupler_comm = MPI_Comm_f2c(comms_shell);
		int cycle_counter = 0;

		int my_rank;
  		MPI_Comm_rank(coupler_comm, &my_rank);

		bool found = false;
		int unit_count = 0;

    	while(!found){
			if(units[unit_count].coupler_ranks.size() > 0) {
				for(int j=units[unit_count].coupler_ranks[0][0]; j<units[unit_count].coupler_ranks[0][0]+units[unit_count].coupler_ranks[0].size();j++){//if rank of processes matches a rank for a particular coupler unit, assign the new communicator
					if((units[unit_count].type == 'C' && rank == j) || (units[unit_count].type == 'L' && rank == j)){
						found=true;
						printf("Found!\n");
					}
				}
			}
            if(found == false){
                unit_count++;
            }
    	}
		
		int left_rank = units[unit_count].mgcfd_ranks[0][0];
		//int left_size = static_cast<int>(units[unit_count].mgcfd_ranks[0].size());
		int right_rank;
		if(units[unit_count].type == 'C'){
			right_rank = units[unit_count].mgcfd_ranks[1][0];
		}

		//int right_size = static_cast<int>(units[unit_count].mgcfd_ranks[1].size());
		//int *left_rank_storage = new int[left_size];
		//int *right_rank_storage = new int[right_size];
        
		int total_ranks = units[unit_count].coupler_ranks[0].size();
		int root_rank = units[unit_count].coupler_ranks[0][0];

		asio::io_context io_context;
        asio::error_code ec;
		asio::ip::tcp::acceptor *acc;
   	    asio::ip::tcp::socket *sock;

		printf("%c\n", units[unit_count].type);

        double left_nodes_sizes[4];
        double right_nodes_sizes[4];
        MPI_Recv(left_nodes_sizes, 4, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(units[unit_count].type == 'C'){
			MPI_Recv(right_nodes_sizes, 4, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else if(units[unit_count].type == 'L'){
			if(rank == root_rank){

				acc = new asio::ip::tcp::acceptor(io_context, asio::ip::tcp::endpoint(asio::ip::tcp::v4(), 65000));
        		sock = new asio::ip::tcp::socket(io_context);

				printf("Listening for a connection\n");
				acc->accept(*sock);
        		printf("Connected!\n");

				std::size_t n = asio::read(*sock, asio::buffer(right_nodes_sizes, 4 * sizeof(double)));
				printf("Sizes received: %f %f %f %f\n", right_nodes_sizes[0], right_nodes_sizes[1], right_nodes_sizes[2], right_nodes_sizes[3]);
			}
			MPI_Bcast(right_nodes_sizes, 4, MPI_DOUBLE, 0, coupler_comm);

		}
       
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

    	double left_nodes_size_chunks[4];
		//double left_extra_chunk[4];
		double right_nodes_size_chunks[4];
		//double right_extra_chunk[4];
		//double left_add_amount[4];
		//double right_add_amount[4];
		//bool chunk_add = true;

		for(int j = 0; j < 4; j++){
				left_nodes_size_chunks[j] = left_nodes_sizes[j] / total_ranks;
				//left_extra_chunk[j] = left_nodes_sizes[j] - left_nodes_size_chunks[j]*total_ranks;//since dividing integers is a floor function in C++, there may be a little bit left which we add onto the first transfer
				right_nodes_size_chunks[j] = right_nodes_sizes[j] / total_ranks;
				//right_extra_chunk[j] = right_nodes_sizes[j] - right_nodes_size_chunks[j]*total_ranks;
            }

        double *left_p_variables_l0_sg, *left_p_variables_l1_sg, *left_p_variables_l2_sg, *left_p_variables_l3_sg;
        double *right_p_variables_l0_sg, *right_p_variables_l1_sg, *right_p_variables_l2_sg, *right_p_variables_l3_sg;

		double vector_counter;
		int sub_count;
		double vector_counter_max;
		double vector_counter_max_sizes_l[4] = {left_nodes_size_chunks[0],left_nodes_size_chunks[1], left_nodes_size_chunks[2], left_nodes_size_chunks[3]};
		double vector_counter_max_sizes_r[4] = {right_nodes_size_chunks[0],right_nodes_size_chunks[1], right_nodes_size_chunks[2], right_nodes_size_chunks[3]};

		left_p_variables_l0_sg = (double *) malloc((left_nodes_size_chunks[0]) * NVAR * sizeof(double)); //left p_variables storage for scatter/gather
        left_p_variables_l1_sg = (double *) malloc((left_nodes_size_chunks[1]) * NVAR * sizeof(double));
        left_p_variables_l2_sg = (double *) malloc((left_nodes_size_chunks[2]) * NVAR * sizeof(double));
        left_p_variables_l3_sg = (double *) malloc((left_nodes_size_chunks[3]) * NVAR * sizeof(double));

        right_p_variables_l0_sg = (double *) malloc((right_nodes_size_chunks[0]) * NVAR * sizeof(double)); //right p_variables storage for scatter/gather
        right_p_variables_l1_sg = (double *) malloc((right_nodes_size_chunks[1]) * NVAR * sizeof(double));
        right_p_variables_l2_sg = (double *) malloc((right_nodes_size_chunks[2]) * NVAR * sizeof(double));
        right_p_variables_l3_sg = (double *) malloc((right_nodes_size_chunks[3]) * NVAR * sizeof(double));

        double *left_p_variable_pointers[4] = {left_p_variables_l0_sg,left_p_variables_l1_sg,left_p_variables_l2_sg,left_p_variables_l3_sg};
        double *right_p_variable_pointers[4] = {right_p_variables_l0_sg,right_p_variables_l1_sg,right_p_variables_l2_sg,right_p_variables_l3_sg};

		double *left_p_variable_pointers_full[4] = {left_p_variables_l0,left_p_variables_l1,left_p_variables_l2,left_p_variables_l3};
        double *right_p_variable_pointers_full[4] = {right_p_variables_l0,right_p_variables_l1,right_p_variables_l2,right_p_variables_l3};

		std::vector< std::vector<double> > vector_of_state_vars_l0, vector_of_state_vars_l1, vector_of_state_vars_l2, vector_of_state_vars_l3; //p_variable vectors to represent brute force search
		std::vector< std::vector<double> > left_vector_of_state_vars_total[4] = {vector_of_state_vars_l0, vector_of_state_vars_l1, vector_of_state_vars_l2, vector_of_state_vars_l3};
		std::vector< std::vector<double> > right_vector_of_state_vars_total[4] = {vector_of_state_vars_l0, vector_of_state_vars_l1, vector_of_state_vars_l2, vector_of_state_vars_l3};

		double comm_total = 0.0;
		double comm_time = 0.0;
		double recv_total = 0.0;
		double recv_time = 0.0;
		double search_total = 0.0;
		double search_time = 0.0;
		double interp_total = 0.0;
		double interp_time = 0.0;

		double coupler_total = 0.0;
		if(rank == root_rank) {
			coupler_total = MPI_Wtime();
		}

		while((cycle_counter < mgcycles) && (cycle_counter % conversion_factor) == 0){// Change this value to the number of cycles if it is not the default

			int local_size;
			MPI_Comm_size(coupler_comm, &local_size);
			if(rank == root_rank){
				MPI_Recv(left_p_variables_l0, left_nodes_sizes[0] * NVAR, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            MPI_Recv(left_p_variables_l1, left_nodes_sizes[1] * NVAR, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            MPI_Recv(left_p_variables_l2, left_nodes_sizes[2] * NVAR, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            MPI_Recv(left_p_variables_l3, left_nodes_sizes[3] * NVAR, MPI_DOUBLE, left_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    	recv_time = MPI_Wtime();

				if(units[unit_count].type == 'C') {
					MPI_Recv(right_p_variables_l0, right_nodes_sizes[0] * NVAR, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            	MPI_Recv(right_p_variables_l1, right_nodes_sizes[1] * NVAR, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            	MPI_Recv(right_p_variables_l2, right_nodes_sizes[2] * NVAR, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	            	MPI_Recv(right_p_variables_l3, right_nodes_sizes[3] * NVAR, MPI_DOUBLE, right_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				else if(units[unit_count].type == 'L') {
					std::size_t n = asio::read(*sock, asio::buffer(right_p_variables_l0, right_nodes_sizes[0] * NVAR * sizeof(double)));
					n = asio::read(*sock, asio::buffer(right_p_variables_l1, right_nodes_sizes[1] * NVAR * sizeof(double)));
					n = asio::read(*sock, asio::buffer(right_p_variables_l2, right_nodes_sizes[2] * NVAR * sizeof(double)));
					n = asio::read(*sock, asio::buffer(right_p_variables_l3, right_nodes_sizes[3] * NVAR * sizeof(double)));
				} 
				recv_time = MPI_Wtime() - recv_time;
				recv_total += recv_time;
	        }
			MPI_Barrier(coupler_comm);
			MPI_Scatter(left_p_variables_l0, (left_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, left_p_variables_l0_sg, (left_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(left_p_variables_l1, (left_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, left_p_variables_l1_sg, (left_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(left_p_variables_l2, (left_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, left_p_variables_l2_sg, (left_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(left_p_variables_l3, (left_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, left_p_variables_l3_sg, (left_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, 0, coupler_comm);

			MPI_Barrier(coupler_comm);
			MPI_Scatter(right_p_variables_l0, (right_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, right_p_variables_l0_sg, (right_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(right_p_variables_l1, (right_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, right_p_variables_l1_sg, (right_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(right_p_variables_l2, (right_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, right_p_variables_l2_sg, (right_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, 0, coupler_comm);
			MPI_Scatter(right_p_variables_l3, (right_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, right_p_variables_l3_sg, (right_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, 0, coupler_comm);

			
			if(MUM == 0){
				MPI_Bcast(left_p_variables_l0, left_nodes_sizes[0] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(left_p_variables_l1, left_nodes_sizes[1] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(left_p_variables_l2, left_nodes_sizes[2] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(left_p_variables_l3, left_nodes_sizes[3] * NVAR, MPI_DOUBLE, 0, coupler_comm);

				MPI_Bcast(right_p_variables_l0, right_nodes_sizes[0] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(right_p_variables_l1, right_nodes_sizes[1] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(right_p_variables_l2, right_nodes_sizes[2] * NVAR, MPI_DOUBLE, 0, coupler_comm);
				MPI_Bcast(right_p_variables_l3, right_nodes_sizes[3] * NVAR, MPI_DOUBLE, 0, coupler_comm);
			}
			
			if (rank == root_rank) {
				search_time = MPI_Wtime();
			}

			//rendezvous routines start
			if(units[unit_count].coupling_type == 'S' || cycle_counter == 0){
				if((cycle_counter % upd_freq) == 0){
					for(int l = 0; l < 5; l++){
						for(int k = 0; k < 4; k++){
				        	vector_counter = 0;
							if(MUM == 0){
								vector_counter_max = std::min(left_nodes_sizes[k], right_nodes_sizes[k]);//this is size of mesh recieved from broadcast
							}else{
								vector_counter_max = std::min(vector_counter_max_sizes_l[k], vector_counter_max_sizes_r[k]);//this is size of mesh recieved from scatter
							}
							sub_count = 0;
							while(sub_count < total_ranks){
								left_vector_of_state_vars_total[k].clear();
								while(vector_counter < (vector_counter_max/total_ranks)){
									std::vector<double> node_state_vars;
									for(int i = 0; i<NVAR; i++){
										if(MUM == 0){
											node_state_vars.push_back(*(left_p_variable_pointers_full[k] + (static_cast<long long>(vector_counter) * NVAR) + i));//essentially move along left_p_variables in chunks of NVAR
										}else{
											node_state_vars.push_back(*(left_p_variable_pointers[k] + (static_cast<long long>(vector_counter) * NVAR) + i));//essentially move along left_p_variables in chunks of NVAR
										}
									}
									left_vector_of_state_vars_total[k].insert(left_vector_of_state_vars_total[k].begin(), node_state_vars);
									vector_counter++;
								}
								vector_counter = 0;
								sub_count++;
							}
			        	}
			        }
				}

				if((cycle_counter % upd_freq) == 0){
					for(int l = 0; l < 5; l++){
						for(int k = 0; k < 4; k++){
				        	vector_counter = 0;
							if(MUM == 0){
								vector_counter_max = std::min(left_nodes_sizes[k], right_nodes_sizes[k]);//this is size of mesh recieved from broadcast
							}else{
								vector_counter_max = std::min(vector_counter_max_sizes_l[k], vector_counter_max_sizes_r[k]);//this is size of mesh recieved from scatter
							}
							sub_count = 0;
							while(sub_count < total_ranks){
								right_vector_of_state_vars_total[k].clear();
								while(vector_counter < (vector_counter_max/total_ranks)){
									std::vector<double> node_state_vars;
									for(int i = 0; i<NVAR; i++){
										if(MUM == 0){
											node_state_vars.push_back(*(right_p_variable_pointers_full[k] + (static_cast<long long>(vector_counter) * NVAR) + i));//essentially move along right_p_variables in chunks of NVAR
										}else{
											node_state_vars.push_back(*(right_p_variable_pointers[k] + (static_cast<long long>(vector_counter) * NVAR) + i));//essentially move along right_p_variables in chunks of NVAR
										}
									}
									right_vector_of_state_vars_total[k].insert(right_vector_of_state_vars_total[k].begin(), node_state_vars);
									vector_counter++;
								}
								vector_counter = 0;
								sub_count++;
							}
			        	}
			        }
				}
			}
			//rendezvous routines end

			if (rank == root_rank) {
				search_time = MPI_Wtime() - search_time;
				search_total += search_time;
			}

			if (rank == root_rank) {
				interp_time = MPI_Wtime();
			}

			//interpolate routine start
			if((cycle_counter % conversion_factor) == 0){
				for(int k = 0; k < 4; k++){
					std::vector<double> node_state_vars_left;
					std::vector<double> node_state_vars_right;
					std::vector<double> node_state_vars_temp;
					//double valuecheck_left;
					//double valuecheck_right;
					vector_counter = 0;
					sub_count = 0;
                    //NEW
                    vector_counter_max = std::min(vector_counter_max_sizes_l[k], vector_counter_max_sizes_r[k]);//this is size of mesh recieved from scatter
					while(sub_count < (total_ranks*500)){
						while(vector_counter < vector_counter_max/total_ranks){
							node_state_vars_left = left_vector_of_state_vars_total[k].at(vector_counter);
							node_state_vars_right = right_vector_of_state_vars_total[k].at(vector_counter);
							node_state_vars_temp = node_state_vars_right;
							for(int i = 0; i<NVAR; i++){
								node_state_vars_right.at(i) = (node_state_vars_left.at(i) + node_state_vars_right.at(i))/2;
								node_state_vars_left.at(i) = (node_state_vars_left.at(i) + node_state_vars_temp.at(i))/2;
							}
							vector_counter++;
						}
						vector_counter = 0;
						sub_count++;
					}
				}
			}
			//interpolate routine end
			
			if (rank == root_rank) {
				interp_time = MPI_Wtime() - interp_time;
				interp_total += interp_time;
			}

			MPI_Barrier(coupler_comm);
	        MPI_Gather(left_p_variables_l0_sg, (left_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, left_p_variables_l0, (left_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(left_p_variables_l1_sg, (left_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, left_p_variables_l1, (left_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(left_p_variables_l2_sg, (left_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, left_p_variables_l2, (left_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(left_p_variables_l3_sg, (left_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, left_p_variables_l3, (left_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, 0, coupler_comm);

			MPI_Barrier(coupler_comm);
	        MPI_Gather(right_p_variables_l0_sg, (right_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, right_p_variables_l0, (right_nodes_size_chunks[0] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(right_p_variables_l1_sg, (right_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, right_p_variables_l1, (right_nodes_size_chunks[1] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(right_p_variables_l2_sg, (right_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, right_p_variables_l2, (right_nodes_size_chunks[2] * NVAR), MPI_DOUBLE, 0, coupler_comm);
	        MPI_Gather(right_p_variables_l3_sg, (right_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, right_p_variables_l3, (right_nodes_size_chunks[3] * NVAR), MPI_DOUBLE, 0, coupler_comm);
            
			if(rank == root_rank){
				comm_time = MPI_Wtime();
	            if(units[unit_count].type == 'C') {
	            	MPI_Send(right_p_variables_l0, right_nodes_sizes[0] * NVAR, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
	            	MPI_Send(right_p_variables_l1, right_nodes_sizes[1] * NVAR, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
	            	MPI_Send(right_p_variables_l2, right_nodes_sizes[2] * NVAR, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
	            	MPI_Send(right_p_variables_l3, right_nodes_sizes[3] * NVAR, MPI_DOUBLE, right_rank, 0, MPI_COMM_WORLD);
				}
				else if(units[unit_count].type == 'L') {
					std::size_t n = asio::write(*sock, asio::buffer(right_p_variables_l0, right_nodes_sizes[0] * NVAR * sizeof(double)));
					n = asio::write(*sock, asio::buffer(right_p_variables_l1, right_nodes_sizes[1] * NVAR * sizeof(double)));
					n = asio::write(*sock, asio::buffer(right_p_variables_l2, right_nodes_sizes[2] * NVAR * sizeof(double)));
					n = asio::write(*sock, asio::buffer(right_p_variables_l3, right_nodes_sizes[3] * NVAR * sizeof(double)));
				}
	            MPI_Send(left_p_variables_l0, left_nodes_sizes[0] * NVAR, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
	            MPI_Send(left_p_variables_l1, left_nodes_sizes[1] * NVAR, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
	            MPI_Send(left_p_variables_l2, left_nodes_sizes[2] * NVAR, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
	            MPI_Send(left_p_variables_l3, left_nodes_sizes[3] * NVAR, MPI_DOUBLE, left_rank, 0, MPI_COMM_WORLD);
				comm_time = MPI_Wtime() - comm_time;
				comm_total += comm_time;
	        }
            
			cycle_counter = cycle_counter + conversion_factor;

		}

		if (rank == root_rank) {
			coupler_total = MPI_Wtime() - coupler_total;
			printf("Total Time (Coupler): %f\n", coupler_total);
			printf("Send Time (Coupler): %f\n", comm_total);
			printf("Recv Time (Coupler): %f\n", recv_total);
			printf("Search Time (Coupler): %f\n", search_total);
			printf("Interp Time (Coupler): %f\n", interp_total);
		}

		MPI_Barrier(coupler_comm);
		MPI_Finalize();
   		exit(0);
	}
}




