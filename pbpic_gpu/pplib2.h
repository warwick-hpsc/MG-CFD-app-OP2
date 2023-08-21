/* header file for pplib2.c */
#include "mpi.h"

void cppinit2(int *idproc, int *nvp, int argc, char *argv[],
			  MPI_Comm pic_comm);

void cppfndgrp(int locl[], int kstrt, int nvp, int *idev, int *ndev);

void cppexit();

void cppabort();

void cpwtimera(int icntrl, float *time, double *dtime);

void cppsum(float f[], float g[], int nxp);

void cppdsum(double f[], double g[], int nxp);

void cppimax(int f[], int g[], int nxp);

void cpppncguard2l(float scs[], float scr[], int kstrt, int nvp,
                   int nxv);

void cpppnaguard2l(float scs[], float scr[], int kstrt, int nvp,
                   int nxv);

void cppptpose(std::complex<float> sm[], std::complex<float> tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp);

void cppptposen(std::complex<float> sm[], std::complex<float> tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim);

void cacsndrec(std::complex<float> stm[], int idproc, int nsize, int ntag,
               int mode);

void cpppmove2(float sbufr[], float sbufl[], float rbufr[], 
               float rbufl[], int ncll[], int nclr[], int mcll[],
               int mclr[], int kstrt, int nvp, int idimp, int nbmax,
               int mx1);
