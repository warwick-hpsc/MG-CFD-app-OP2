/* header file for pplib2.c */
#include "mpi.h"
#include <stdbool.h>

void cppinit2(int *idproc, int *nvp, int argc, char *argv[],
			  MPI_Comm pic_comm);

void cppexit();

void cppabort();

void cpwtimera(int icntrl, float *time, double *dtime);

void cppsum(float f[], float g[], int nxp);

void cppdsum(double f[], double g[], int nxp);

void cppimax(int f[], int g[], int nxp);

void cppdmax(double f[], double g[], int nxp);

void cppncguard2l(float f[], int nyp, int kstrt, int nvp, int nxv,
                  int nypmx);

void cppnaguard2l(float f[], float scr[], int nyp, int nx, int kstrt,
                  int nvp, int nxv, int nypmx);

void cppnacguard2l(float f[], float scr[], int nyp, int nx, int ndim,
                   int kstrt, int nvp, int nxv, int nypmx);

void cpptpose(std::complex<float> f[], std::complex<float> g[], 
			  std::complex<float> s[], std::complex<float> t[], 
			  int nx, int ny, int kxp, int kyp, int kstrt, 
			  int nvp, int nxv, int nyv, int kxpd, int kypd);

void cppntpose(std::complex<float> f[], std::complex<float> g[], 
			   std::complex<float> s[], std::complex<float> t[], 
			   int nx, int ny, int kxp, int kyp, int kstrt, 
			   int nvp, int ndim, int nxv, int nyv, int kxpd, int kypd);

void cppmove2(float part[], float edges[], int *npp, float sbufr[],
              float sbufl[], float rbufr[], float rbufl[], int ihole[],
              int ny, int kstrt, int nvp, int idimp, int npmax, int idps,
              int nbmax, int ntmax, int info[], FILE *fp, bool debug);
