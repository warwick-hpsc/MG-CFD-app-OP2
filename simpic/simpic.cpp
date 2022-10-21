#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>


#include "global.h"
#include <mpi.h>
#include "timing.cpp"
#include "diagnostics.cpp"
#include "trimatrix.cpp"
#include "fields.cpp"
#include "coupler_config.h"
#include "../src/structures.h"


inline void lhs(int j)
{
  lhsbuf[nlp] = pdata[j];
  lhsbuf[nlp+1] = pdata[j+1];
  nlp += 2;

#ifdef DEBUG
  fprintf(stderr, "rank %d nlp => %d\n", rank, nlp);
#endif 

}

inline void rhs(int j)
{
  rhsbuf[nrp] = pdata[j];
  rhsbuf[nrp+1] = pdata[j+1];
  nrp += 2;

#ifdef DEBUG
  fprintf(stderr, "rank %d nrp => %d\n", rank, nrp);
#endif 
}

inline Scalar gridx(Scalar x)
{
  return ((x - xl)/dx);
}

inline Scalar getfrac(Scalar x)
{
  static int n;

  n = int(x);
  return x-n;
}

inline void weightparticle(Scalar x, Scalar weight, Scalar *array)
{
  
  xx = gridx(x);
  n = int(xx);
  frac = xx-n;

  array[n] += weight*(1. - frac);
  array[n+1] += weight*frac;
}

inline Scalar particleweight(Scalar x, Scalar *array)
{

  xx = gridx(x);
  n = int(xx);
  frac = xx-n;

  return frac*array[n+1] + (1.0 - frac)*array[n];
}

void moveparticles(Scalar dt)
{

  starttime(PARTICLES);
  for(i=0, j=0; i < npart; )
    {
      E = particleweight(pdata[j], Earray);
      pdata[j+1] += qm*E;
      pdata[j] += pdata[j+1]*dt;

      #ifdef DEBUG
      //      fprintf(stderr, "particle %d (%d): E=%g, x=%g, v=%g, qm=%g\n", i, j,
      //	      E, pdata[j], pdata[j+1], qm);
      #endif

      if((pdata[j] > xl) && (pdata[j] < xr))
	{
	  weightparticle(pdata[j], qscale, narray);
	  j += 2;
	  i++;
	}
      else
	{
	  #ifdef DEBUG
	  //  fprintf(stderr, "\nParticle %d[of %d] lost, E=%g\n", i, npart, E);
	  #endif 

	  if(pdata[j] < xl)
	    {
	      lhs(j);
	    }
	  else
	    {
	      rhs(j);
	    }	  
	  npart--;
	  pdata[j] = pdata[2*npart];
	  pdata[j+1] = pdata[2*npart + 1];
	}
    }
  endtime(PARTICLES);
}

void communicateparticles(void)
{

#ifdef DEBUG
  fprintf(stderr, 
	  "Rank %d communicateparticles(), t=%g, nlp=%d, nrp=%d, npart=%d\n", 
	  rank, tt, nlp, nrp, npart);
  
#endif

  // recv/send from/to left
  if(rank != 0)
    {

#ifdef DEBUG
      fprintf(stderr, "left: Processor %d sending %d particles to proc %d\n",
	      rank, nlp, nl);
#endif

      MPI_Irecv(&nlr, 1, MPI_INT, nl, 0, custom_comm, &lstat); 
      MPI_Isend(&nlp, 1, MPI_INT, nl, 0, custom_comm, &llstat);
      MPI_Request_free(&llstat);
    }

  // recv/send from/to right
  if(rank != last)
    {

#ifdef DEBUG
      fprintf(stderr, "right: Processor %d sending %d particles to proc %d\n",
	      rank, nrp, nr);
#endif

      MPI_Irecv(&nrr, 1, MPI_INT, nr, 0, custom_comm, &rstat);
      MPI_Isend(&nrp, 1, MPI_INT, nr, 0, custom_comm, &rrstat);
      MPI_Request_free(&rrstat);
    }
      
#ifdef DEBUG
  fprintf(stderr, "proc %d waiting for # of particles, t=%g\n", rank, tt);
#endif

  if(rank != 0)
    {
#ifdef DEBUG
      fprintf(stderr, "proc %d: rank != 0, t=%g\n", rank, tt);
#endif

      MPI_Wait(&lstat, &stat);
      //MPI_Request_free(&lstat);  
#ifdef DEBUG
      fprintf(stderr, "left: Processor %d receiving %d particles from proc %d\n",
	      rank, nlr, nl);
#endif
      //MPI_Request_free(&lstat);
      //MPI_Request_free(&llstat);
    }
  if(rank != last)
    {
#ifdef DEBUG
      fprintf(stderr, "proc %d: rank != last, t=%g\n", rank, tt);
#endif

      MPI_Wait(&rstat, &stat);
#ifdef DEBUG
      fprintf(stderr, "right: Processor %d receving %d particles from proc %d\n",
	      rank, nrr, nr);
#endif
      //MPI_Request_free(&rstat);
      //MPI_Request_free(&rrstat);
    }
  
#ifdef DEBUG
  fprintf(stderr, "Rank %d ready for particle communication, t=%g:\n", rank, tt);
#endif


  // here, one would reallocate, if necessary
  // for this simple test, we use unrealistically
  // large buffers to avoid this possibility

  if(rank != 0)
    {
      // receive from left
      MPI_Irecv(&(pdata[2*npart]), nlr, MPI_SCALAR, nl, 0, custom_comm, &plstat);

      npart += (nlr/2);

      // send to left
      MPI_Isend(lhsbuf, nlp, MPI_SCALAR, nl, 0, custom_comm, &pllstat);
      MPI_Request_free(&plstat);
      MPI_Request_free(&pllstat);

      //      if(rank != last)
      //	{
      //      MPI_Wait(&plstat, &stat);
	  //	}
      //MPI_Request_free(&plstat);
      //MPI_Request_free(&pllstat);
    }
  if(rank != last)
    {
      MPI_Irecv(&(pdata[2*npart]), nrr, MPI_SCALAR, nr, 0, custom_comm, &prstat);
      npart += (nrr/2);

      // send to right
      MPI_Isend(rhsbuf, nrp, MPI_SCALAR, nr, 0, custom_comm, &prrstat);
      MPI_Request_free(&prrstat);

      MPI_Wait(&prstat, &stat);
      //MPI_Request_free(&prstat);
      //MPI_Request_free(&prrstat);
    }
  //  if(rank != 0) { MPI_Wait(&plstat, &stat); }

  #ifdef DEBUG
  fprintf(stderr, "rank %d waiting for MPI_Barrier(), t=%g\n", rank, tt);
  #endif

  // wait for all processors to catch up
  MPI_Barrier(custom_comm);

  #ifdef DEBUG
  fprintf(stderr, "rank %d done with communicateparticles(), t=%g.\n", rank, tt);
  #endif

}

void print_diagnostics(char *fmod, int close=0)
{
  static int init = 0;
  char ffname[STR_SIZE];
  char temp[STR_SIZE];

  if(!init)
    {
      #ifdef DEBUG
      //      fprintf(stderr, "\nInitializing diagnostics...\n");
      #endif

      sprintf(diags[NDENSITY], "density");
      thediags[NDENSITY] = narray;
      sprintf(diags[NEFIELD], "E");
      thediags[NEFIELD] = Earray;
      sprintf(diags[NVOLTAGE], "phi");
      thediags[NVOLTAGE] = phiarray;
      sprintf(diags[NVXPHASE], "vxx");
      thediags[NVXPHASE] = pdata;

      for(int i=0; i < NDIAGNOSTICS; i++)
	{
	  sprintf(temp, "%s%s", diags[i], procname);
	  strcpy(diags[i], temp);
	}

      init = 1;
    }
  
  if(close)
    {
      return;
    }

  for(int i=0; i < NDIAGNOSTICS; i++)
    {
      sprintf(ffname, "%s%s%s", diags[i], fmod, fext);
      if(i < NVXPHASE)
	{
	  printtarray(i, ng, dx);
	}
      else
	{
	  printphase(ffname);	  
	}
    }
  printt();
}

void mainloop(Scalar tmax, Scalar dt)
{
  tt = 0.;
  int count = 0;
  int comm_size = 0;
  double interface_percentage = 0.05;
  MPI_Comm_size(custom_comm, &comm_size);
  Scalar *narray_variables = NULL;
  narray_variables = new Scalar[ng * comm_size];
  Scalar *large_interface;
  Scalar *large_interface_recv;
  Scalar transfer_size;
  double interface_size;
  int coupler_rank;

  if(rank == 0){
    interface_size = std::round(0.05 * 1000000 * artificalsize);
    large_interface = new Scalar[interface_size];
    large_interface_recv = new Scalar[interface_size];
    std::fill_n(large_interface, interface_size, 0);
    if(ng * comm_size > interface_size){
      transfer_size = interface_size;
    }else{
      transfer_size = ng * comm_size;
    }
    int ranks_per_coupler;
    for(int z = 0; z < total_coupler_unit_count; z++){
      ranks_per_coupler = units_copy[unit_count].coupler_ranks[z].size();
      for(int z2 = 0; z2 < ranks_per_coupler; z2++){
          MPI_Send(&interface_size, 1, MPI_DOUBLE, units_copy[unit_count].coupler_ranks[z][z2], 0, MPI_COMM_WORLD);//this sends the node sizes to each of the coupler ranks of each of the coupler units
      }
    }
  }
  
  while(tt < tmax)
    {
      if(count % (ntimesteps/coupler_cycles) == 0 && count < (ntimesteps/coupler_cycles) * coupler_cycles){//this will ensure coupling takes place the right number of times - note that coupling doesn't take place on the final iteration
        MPI_Barrier(custom_comm);
        MPI_Gather(narray, ng, MPI_SCALAR, narray_variables, ng, MPI_SCALAR, 0, custom_comm);//gather the SIMPIC mesh data from each of the ranks
        if(rank == 0){
          printf("Count is %d, sending from simpic side\n", count);  
          std::memcpy(large_interface, narray, transfer_size);
          for(int z = 0; z < total_coupler_unit_count; z++){
            coupler_rank = units_copy[unit_count].coupler_ranks[z][0];
            MPI_Send(large_interface, interface_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD);
            MPI_Recv(large_interface_recv, interface_size, MPI_DOUBLE, coupler_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          printf("Count is %d, receiving from simpic side\n", count); 
        }
        MPI_Barrier(custom_comm);
      }
      
      #ifdef DEBUG
      fprintf(stderr, "rank %d mainloop: t=%g, npart=%d{\n", 
	      rank, tt, npart);
      #endif

      // reset needed variables
      nlp = nrp = 0;
      memcpy(narray, nback, ng*sizeof(Scalar));

      nparttot += npart;

      moveparticles(dt);
      
      #ifdef DEBUG
      //      fprintf(stderr, "after move particles, npart=%d\n", npart);
      #endif

      if(nproc > 1)
	{
	  starttime(MPITIME);
	  communicateparticles();
	  endtime(MPITIME);
	}
      

      advancefields(dt);

      #ifdef DEBUG
      //      fprintf(stderr, "after fields, npart=%d\n", npart);
      #endif

      if(diagnosticsflag)
	{
	  char msg[] = "";
	  starttime(DIAGNOSTICS);
	  print_diagnostics(msg);
	  endtime(DIAGNOSTICS);
	}

      #ifdef DEBUG
      //      fprintf(stderr, "after diagnostics, npart=%d\n", npart);
      #endif

      tt += dt;

      MPI_Barrier(custom_comm);

      #ifdef DEBUG
      fprintf(stderr, "}[%d]\n", rank);
      #endif
      count++;
    }
    delete [] narray_variables;
}

void allocate_particles(void)
{
  int j;
  int i;
  // nproc + 1 -- safety factor for this test
  // (should be more dynamic)
  pdata = new Scalar[(10*npart+1)];
  lhsbuf = new Scalar[(10*npart+1)];
  rhsbuf = new Scalar[(10*npart+1)];

  for(i=0, j=0; i < npart; i++, j+=2)
    {
      pdata[j] = xl + ((i + 1)*(xr-xl))/Scalar(npart + 1);
      pdata[j+1] = 0.;
    }
}



void free_particles(void)
{
  delete [] pdata;
  delete [] lhsbuf;
  delete [] rhsbuf;
}

void init(void)
{

  // remove old .dat files
  system("rm -f *.sdat");

  // get MPI data
  MPI_Get_version(&ver, &subver);  
  MPI_Comm_size (custom_comm, &nproc);
  MPI_Comm_rank (custom_comm, &rank);
  wtres = MPI_Wtick();

  // print processor string name
  sprintf(procname, "%03d", rank);

  // reset timer variables
  for(int j=0; j < NTIMES; j++)
    {
      tttimes[j] = 0.;
    }
  // label timer counters -- do this here because
  // I can't find any documentation on how to declare 
  // multidimensional variables at compile time
  // aka C is stupid
  sprintf(names[PARTICLES], "Particles");
  sprintf(names[FIELDS], "Fields");
  sprintf(names[MPITIME], "MPI");
  sprintf(names[DIAGNOSTICS], "Diagnostics");
  sprintf(names[INIT], "Initialization"); 

  nparttot = 0.;

  //  ntimesteps = 4; // now obtained from command line
  density = 1.E13;
  //  np2c = 1E7; // now obtained from command line
  epsilon = 8.85E-12;
  area = 1.;
  L = 1.;
  q = 1.602E-19;
  m = 9.1E-31;
  //  ng = NG; // now obtained from command line
  //  nc = ng - 1; // now obtained from command line

  xl = rank/Scalar(nproc);
  xr = (rank + 1.)/Scalar(nproc);
  xl *= L;
  xr *= L;
  
  Ll = xr - xl;
  dx = Ll / Scalar(nc);

  nl = rank - 1;
  nr = rank + 1;
  last = nproc - 1;
  if(rank == 0)
    {
      nl = last;
    }
  if(rank == last)
    {
      nr = 0;
    }

  //  npart = int(density*area*Ll/np2c);
  np2c = density*area*Ll/Scalar(npart);

  allocate_particles();
  init_fields();

  // calculate time step from plasma frequency
  wp = density*q*q/(epsilon*m);
  wp = sqrt(wp);

  dt = dtfactor/wp;

  qscale = np2c/(area*dx);

  t = ntimesteps*dt;
  qm = q*dt/m;

  diagnosticsflag = false;
}

void printdata(void)
{
  printf("rank %d of %d processors\n", rank, nproc);
  fprintf(stdout, "MPI version: %d.%d, ", ver, subver);
  fprintf(stdout, "MPI_Wtime precision: %g\n", wtres);
  printf("t=%.3g, dt=%.3g, ntimesteps=%d, ", t, dt, ntimesteps);
  printf("density=%.3g, wp=%.3g, np2c=%.3g\n", density, wp, np2c);
  printf("%d particles, ng=%d, nc=%d, ", npart, ng, nc);
  printf("L=%g, Ll=%g, dx=%g, area=%g\n", L, Ll, dx, area);

  fflush(stdout);
}

// free fields arrays
void free_fields(void)
{
  allocate_arrays(1);
}

void parsefromfile()
{
    int i;
    int nnppcc;   // number of particles per cell
    int nnccpppp; // number of cells per processor
    int nnnttt;   // number of time steps
    int asz; // artificial mesh size for coupling
    FILE *fp;
    char filebuffer[100]; // buffer to hold contents from simpic input file
    int count; // used to count the total number of input parameters (i.e to replace argc)
    double ddttff, lhsdefault;
    int array_count; // used to store the current array index

    count = array_count = 0;
    ddttff = lhsdefault = 0.;
    nnppcc = nnccpppp = nnnttt = 0;

    fp = fopen("simpic_parameters", "r");
    if(fp == NULL){
        printf("Error: SIMPIC parameter file not found\n");
        Quit();
    }

    if(filebuffer != fgets(filebuffer, 100, fp)){
        printf("SIMPIC parameter file exists but something went wrong\n");
        Quit();
    }

    for (i = 0; filebuffer[i] != '\0'; i++){
        if (filebuffer[i] == ' ' && filebuffer[i+1] != ' '){
            count++;
        }
    }
    count++; //increment count since the last word will not have a space after


    char* word = strtok(filebuffer," ");
    char** fiparameters = (char**) calloc(count, sizeof(char*)); //allocate an array of strings
    fiparameters[array_count] = (char*) calloc(strlen(word) + 1, sizeof(char)); //allocate memory to an individual string
    strncpy(fiparameters[array_count], word, strlen(word)+1);
    array_count++;
    while(word != NULL) {
        word = strtok(NULL, " ");
        if (word != NULL) {
            fiparameters[array_count] = (char*) calloc(strlen(word) + 1, sizeof(char)); //allocate memory to the next individual string
            strncpy(fiparameters[array_count], word, strlen(word)+1);
            array_count++;
        }
    }

    for(i=0; i < count; i++)
    {

        if(strcmp(fiparameters[i], "-ppc") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%d", &nnppcc);
        }
        else if(strcmp(fiparameters[i], "-asz") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%d", &asz);
        }
        else if(strcmp(fiparameters[i], "-ncpp") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%d", &nnccpppp);
        }
        else if(strcmp(fiparameters[i], "-nt") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%d", &nnnttt);
        }
        else if(strcmp(fiparameters[i], "-dtfactor") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%lf", &ddttff);
        }
        else if(strcmp(fiparameters[i], "-lhsv") == 0)
        {
            i++;
            sscanf(fiparameters[i], "%lf", &lhsdefault);
        }
        else // default case
        {
            fprintf(stderr, "\nError:\n");
            fprintf(stderr, "Unrecognized argument \"%s\"\n", fiparameters[i]);
            Quit();
        }
    }


    if((nnppcc < 1) || (nnccpppp < 1) || (nnnttt < 1))
    {
        fprintf(stderr, "\nError, input arguments must be entered!\n");
        Quit();
    }

    if(ddttff <= 0.)
    {
        fprintf(stderr, "\nError, dtfactor MUST BE positive!\n");
        Quit();
    }

    // set simulation variables from input data
    artificalsize = asz;
    ntimesteps = nnnttt;
    nc = nnccpppp;
    ng = nc + 1;
    npart = nc*nnppcc;
    dtfactor = ddttff;
    lhsvoltage = lhsdefault;

    /*
    printf("Number of timesteps: %d\n", ntimesteps);
    printf("Number of cells per processor: %d\n", nc);
    printf("Number of particles per cell: %d\n", nnppcc);
    printf("Dtfactor: %f\n", dtfactor);
    printf("Left hand side voltage: %f\n", lhsvoltage);
    */

    for(i=0; i < count; i++){//free each of the parameter sub-arrays
        free(fiparameters[i]);
    }
    free(fiparameters);//free the parameter array

    fclose(fp);
}

void parsecmdline(int argc, char **argv)
{
  int i;
  int nnppcc;   // number of particles per cell
  int nnccpppp; // number of cells per processor 
  int nnnttt;   // number of time steps
  double asz; // artificial mesh size for coupling
  double ddttff, lhsdefault;

  #ifdef DEBUG
  //  fprintf(stderr, "Program name: %s\n", argv[0]);
  #endif

  ddttff = lhsdefault = 0.;
  nnppcc = nnccpppp = nnnttt = 0;

  for(i=1; i < argc; i++)
    {
      #ifdef DEBUG
      //      fprintf(stderr, "arg %d : %s\n", i, argv[i]);
      #endif

      if(strcmp(argv[i], "-ppc") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnppcc);
	}
      else if(strcmp(argv[i], "-ncpp") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnccpppp);
	}
      else if(strcmp(argv[i], "-asz") == 0)
  {
    i++;
    sscanf(argv[i], "%lf", &asz);
  }
      else if(strcmp(argv[i], "-nt") == 0)
	{
	  i++;
	  sscanf(argv[i], "%d", &nnnttt);
	}
      else if(strcmp(argv[i], "-dtfactor") == 0)
	{
	  i++;
	  sscanf(argv[i], "%lf", &ddttff);
	}
      else if(strcmp(argv[i], "-lhsv") == 0)
	{
	  i++;
	  sscanf(argv[i], "%lf", &lhsdefault);
	}
      else // default case
	{
	  fprintf(stderr, "\nError:\n");
	  fprintf(stderr, "Unrecognized argument \"%s\"\n", argv[i]);
	  Quit();
	}
    }

  #ifdef DEBUG
  //  fprintf(stderr, "\ncmd line args parsed: ");
  //  fprintf(stderr, "%d particles/cell, %d cells/processor, %d timesteps\n",
  //	  nnppcc, nnccpppp, nnnttt);
  #endif

  if((nnppcc < 1) || (nnccpppp < 1) || (nnnttt < 1))
    {
      fprintf(stderr, "\nError, input arguments must be entered!\n");
      Quit();
    }

  if(ddttff <= 0.)
    {
      fprintf(stderr, "\nError, dtfactor MUST BE positive!\n");
      Quit();
    }

  // set simulation variables from input data
  artificalsize = asz;
  ntimesteps = nnnttt;
  nc = nnccpppp;
  ng = nc + 1;
  npart = nc*nnppcc;
  dtfactor = ddttff;
  lhsvoltage = lhsdefault;
}

void Quit(void)
{
  free_particles();
  free_fields();

  allocate_fieldbuffers(1);

  if(diagnosticsflag)
    {
      printt(1);
      printtarray(0, 0, 0., 1);
    }

  #ifdef DEBUG
  //  fprintf(stderr, "proc %d ready to finalize MPI\n", rank);
  #endif


  //exit(0);
  //MPI_Abort(custom_comm, 0);

   MPI_Finalize();
   exit(0);
}



int main_simpic(int argc, char** argv, MPI_Fint custom, int instance_number, struct unit units[], struct locators relative_positions[])
{
  //MPI_Init(&argc, &argv);
  int flag = 0;
  MPI_Initialized(&flag);
  if (!flag) {
    MPI_Init(&argc, &argv);
  }

  custom_comm = MPI_Comm_f2c(custom);
  //custom_comm = MPI_COMM_WORLD;
  tstart = MPI_Wtime();

  starttime(INIT);

  //parsecmdline(argc, argv);

  parsefromfile();
  init();
  
  #ifdef DEBUG
  if(rank == 0)
    {
      printdata();
    }
  #endif

  endtime(INIT);

  /* Here is where the coupler rank is worked out*/
  int worldrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
  int simpic_unit_num = relative_positions[worldrank].placelocator;
  unit_count = 0;
  int simpic_count = 1;//since units start from 1
  bool found = false;
  while(!found){
      if(units[unit_count].type == 'P' && simpic_unit_num == simpic_count){
          found=true;
      }else {
          if(units[unit_count].type != 'C'){
              simpic_count++;
          }
          unit_count++;
      }
  }
  total_coupler_unit_count = units[unit_count].coupler_ranks.size();
  int coupler_rank = units[unit_count].coupler_ranks[0][0]; //This assumes only 1 coupler unit per 2 MG-CFD sessions
  if(rank == 0){
    printf("Number of ranks in my coupler unit: %zu \n", units[unit_count].coupler_ranks[0].size());
    printf("Number of SIMPIC iterations: %d\n", ntimesteps);
    printf("Number of coupler cycles: %d\n", coupler_cycles);
    printf("Communication happening every %d SIMPIC iterations\n", ntimesteps/coupler_cycles);
    printf("Artificial interface size: %lf million", artificalsize);
  }

  units_copy = units;
  relative_positions_copy = relative_positions;


  mainloop(t, dt);

  #ifdef DEBUG
  //  fprintf(stderr, "Proc %d done with mainloop\n", rank);
  #endif
  MPI_Barrier(custom_comm);

  tend = MPI_Wtime();
  
  if(rank == 0)
    {
      printtimes(stdout);
    }

  #ifdef DEBUG
  //  fprintf(stderr, "Proc %d ready to Quit()\n", rank);
  #endif

  Quit();

#ifdef DEBUG
  //  fprintf(stderr, "Proc %d DONE!!!\n", rank);
#endif
  return 0;
}
