int upd_freq = 810;/* in every x MG cycles, send data through CPX*/
int mgcycles = 800;/* number of MG cycles */
int conversion_factor = 10; //approx. how many MG cycles equals one production cycle - used to control interpolation routine
int MUM = 0; //multi-unit mode - 0 is single unit, multi rank, 1 is multi unit
