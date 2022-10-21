// global variables

#include <mpi.h>
#include "../src/structures.h"
#define Scalar double
#define MPI_SCALAR MPI_DOUBLE

#define STR_SIZE 501

// #define DTFACTOR 0.05  -- no longer used

// #define DEBUG

#define NG 100
#define NP 10

// definitions for diagnostics
#define NDENSITY 0
#define NEFIELD  1
#define NVOLTAGE 2
#define NVXPHASE 3
#define NDIAGNOSTICS 4

Scalar dtfactor;

Scalar *pdata;
Scalar *lhsbuf, *rhsbuf;
int npart, lpart, rpart;

MPI_Comm custom_comm;
// fields arrays
Scalar *Earray;
Scalar *phiarray;
Scalar *narray;
Scalar *nback;

Scalar area;
Scalar density;
Scalar np2c;
Scalar q, m, qm;
Scalar qscale;
Scalar epsilon;
Scalar wp;
Scalar El, Er;
Scalar dx, dt, t;
Scalar xl, xr, L, Ll;
int ntimesteps;
Scalar artificalsize;
struct unit *units_copy;
struct locators *relative_positions_copy;

int total_coupler_unit_count;
int unit_count;
int lproc, rproc;
int nproc, rank;
int ver, subver;
int last;

int ng, nc, ngglobal;

int nl, nr;
int nlp, nrp;  // counters for sending particle buffers
int nlr, nrr;

bool diagnosticsflag;

Scalar lhsvoltage;

// variables for timing
#define PARTICLES 0
#define FIELDS 1
#define MPITIME 2
#define DIAGNOSTICS 3
#define INIT 4
#define NTIMES 5
Scalar tttimes[NTIMES];
char names[NTIMES][STR_SIZE];
Scalar tstart, tend;  // global times
Scalar wtres;
Scalar nparttot;

Scalar tt;  // universal time

// variables for field communication
Scalar *frbuffer;
Scalar *fsbuffer;

char diags[NDIAGNOSTICS][STR_SIZE];
Scalar * thediags[NDIAGNOSTICS];

char fext[] = ".sdat";

char procname[STR_SIZE];

// communicateparticles variables
MPI_Request lstat, rstat;
MPI_Request plstat, prstat;
MPI_Request llstat, rrstat;
MPI_Request pllstat, prrstat;
MPI_Status stat;

// moveparticles variables
int i, j;
Scalar E, x, v;

// weightparticle and particleweight variables
int n;
Scalar frac, xx;

// global functions
void Quit(void);

inline Scalar xcoord(int i)
{
  return xl + i*dx;
}

void gradient(Scalar *grad, Scalar *arr, int n, Scalar scale)
{
  register int i;

  n--;
  for(i=1; i < n; i++)
    {
      grad[i] = scale*(arr[i+1] - arr[i-1]);
    }
}
