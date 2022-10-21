static int coupler_cycles = 5000;/* number of communication cycles */
static int mg_conversion_factor = 10; //approx. how many MG cycles equals one production cycle - how many MG cycles run for each coupler cycle.
static int fenics_conversion_factor = 1; //approx. how many FENICS cycles equals one production cycle - how many FENICS cycles run for each coupler cycle. 
static int search_freq = 6;/* in every x coupler cycles, run the search routine in cpx*/
static int MUM = 1; //multi-unit mode - 0 is single unit, multi rank, 1 is multi unit
static bool fastsearch = true; //when true, tree based search is used, else brute force search is used
static bool ultrafastsearch = true; //mimics the effects of a 'next cell' prediction feature 
static bool superdebug = false; // Disables coupling entirely and allows applications to run on their own
static bool debug = false; //controls the amount of output from cpx
static bool hide_search = false; //controls the overlapping of MG-CFD and cpx Search..
