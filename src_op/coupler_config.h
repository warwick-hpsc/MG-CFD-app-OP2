int upd_freq = 10;/* in every x MG cycles, send data through CPX*/
int mgcycles = 50;/* number of MG cycles */
int conversion_factor = 10; //approx. how many MG cycles equals one production cycle - used to control interpolation routine
int MUM = 0; //multi-unit mode - 0 is single unit, multi rank, 1 is multi unit
bool fastsearch = false; //controls the search type in CPX
bool superdebug = false; // Disables coupling entirely and allows applications to run on their own
bool debug = true; //controls the amount of output from cpx
bool hide_search = false; //controls the overlapping of MG-CFD and cpx Search..
