int upd_freq = 60;/* in every x MG cycles, send data through CPX*/
int mgcycles = 250;/* number of MG cycles */
int conversion_factor = 1; //approx. how many MG cycles equals one production cycle - used to control interpolation routine
int MUM = 1; //multi-unit mode - 0 is single unit, multi rank, 1 is multi unit
