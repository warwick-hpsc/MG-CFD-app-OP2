void send_num_data(struct unit units[], struct locators relative_positions[],
                   long long nodes_size, double **p_variables_data, 
				   double **p_variables_recv);

double send_recv_data(struct unit units[], struct locators relative_positions[],
                    long long nodes_size, int cycle_num, int total_cycles, 
					double *p_variables_data, double *p_variables_recv);

void read_inputs(int *indx, int *indy, int *npx, int *npy);
