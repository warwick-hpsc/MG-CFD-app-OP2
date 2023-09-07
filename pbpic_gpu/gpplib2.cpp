/*--------------------------------------------------------------------*/
/* Basic parallel PIC library for GPU-MPI communications
   gpplib2.c contains basic communications procedures for 1d partitions:
   cgppcacguard2l accumulates guard cells and copies to 3 component vector
                  field
   cgppcaguard2l accumulates guard cells and copies to scalar field
   cgppcbguard2l replicates guard cells and copies to 3 component vector
                field
   cgpporder2l sorts partiles by tiles
   cgpptpose performs a transpose of a complex scalar array, distributed
             in y, to a complex scalar array, distributed in x.
             data from GPU is sent asynchronous, overlapping with MPI
   cgpptposen performs a transpose of an n component complex vector array,
              distributed in y, to an n component complex vector array,
              distributed in x.
              data from GPU is sent asynchronous, overlapping with MPI
   written by viktor k. decyk, ucla
   copyright 2013, regents of the university of california
   update: april 28, 2014                                              */

/*#include <complex>*/
#include <cuda/std/complex>
#include <sys/time.h>
#include "gpulib2.h"
#include "gpupbpush2.h"
#include "gpupfft2.h"
#include "pplib2.h"
#include "gpplib2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

/*--------------------------------------------------------------------*/
void cgppcacguard2l(cuda::std::complex<float> g_cu[], float g_cue[], float g_scs[],
                    float scs[], float scr[], int nx, int nyp,
                    int kstrt, int nvp, int ndim, int nxe, int nypmx,
                    int nxvh, int kypd) {
/* this subroutine copies vector 3 component field and accumulates */
/* guard cells in y from remote GPU into vector field              */
   cgpuppcacguard2xl(g_cu,g_scs,g_cue,nyp,nx,nxe,nypmx,nxvh,kypd);
   gpu_fcopyout(scs,g_scs,2*ndim*nxvh);
   cpppnaguard2l(scs,scr,kstrt,nvp,2*ndim*nxvh);
   gpu_fcopyin(scr,g_scs,2*ndim*nxvh);
   cgpuppcacguard2yl(g_cu,g_scs,nx,nxvh,kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cgppcaguard2l(cuda::std::complex<float> g_q[], float g_qe[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int nxe, int nypmx, int nxvh, int kypd) {
/* this subroutine copies scalar field and accumulates guard cells */
/* in y from remote GPU into scalar field                          */
   cgpuppcaguard2xl(g_q,g_scs,g_qe,nyp,nx,nxe,nypmx,nxvh,kypd);
   gpu_fcopyout(scs,g_scs,2*nxvh);
   cpppnaguard2l(scs,scr,kstrt,nvp,2*nxvh);
   gpu_fcopyin(scr,g_scs,2*nxvh);
   cgpuppcaguard2yl(g_q,g_scs,nx,nxvh,kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cgppcbguard2l(cuda::std::complex<float> g_fxyz[], float g_fxyze[],
                   float g_scs[], float scs[], float scr[], int nx, 
                   int nyp, int kstrt, int nvp, int ndim, int nxe,
                   int nypmx, int nxvh, int kypd) {
/* this subroutine copies 3 component vector field and adds additional */
/* guard cells in y from remote GPU into extended vector field         */
/* local data */
   int nnxe;
   nnxe = ndim*nxe;
   cgpuppcbguard2xl(g_fxyz,g_scs,g_fxyze,nyp,nx,nxe,nypmx,nxvh,kypd);
   gpu_fcopyout(scs,g_scs,2*nxvh*ndim);
   cpppncguard2l(scs,scr,kstrt,nvp,nnxe);
   gpu_fcopyin(scr,g_scs,2*nxvh*ndim);
   cgpuppcbguard2yl(g_fxyze,g_scs,nyp,nx,nxe,nxvh,nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cgpporder2l(float g_ppart[], float g_ppbuff[], float g_sbufl[],
                 float g_sbufr[], int g_kpic[], int g_ncl[],
                 int g_ihole[], int g_ncll[], int g_nclr[],
                 float sbufl[], float sbufr[], float rbufl[],
                 float rbufr[], int ncll[], int nclr[], int mcll[],
                 int mclr[], float ttp[], int noff, int nyp, int kstrt,
                 int nvp, int idimp, int nppmx, int nx, int ny, int mx,
                 int my, int mx1, int myp1, int npbmx, int ntmax,
                 int nbmax, int *g_irc) {
/* this subroutine performs an mpi-gpu particle sort by x,y grid in tiles
   of mx, my
   linear interpolation, with periodic boundary conditions
   for distributed data, with 1d domain decomposition in y.
local data */
   float time;
   double dtime;
   struct timeval itime;
/* first part of particle reorder on x and y cell with mx, my tiles */
   dtimer(&dtime,&itime,-1);
   cgpupppord2la(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,g_ihole,
                 g_ncll,g_nclr,noff,nyp,idimp,nppmx,nx,ny,mx,my,mx1,
                 myp1,npbmx,ntmax,nbmax,g_irc);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[0] += time;
/* move particles on GPU into appropriate spatial regions */
   dtimer(&dtime,&itime,-1);
   gpu_icopyout(ncll,g_ncll,3*mx1);
   gpu_icopyout(nclr,g_nclr,3*mx1);
   gpu_fcopyout(sbufl,g_sbufl,idimp*ncll[3*mx1-1]);
   gpu_fcopyout(sbufr,g_sbufr,idimp*nclr[3*mx1-1]);
   cpppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,
             idimp,nbmax,mx1);
   gpu_icopyin(mcll,g_ncll,3*mx1);
   gpu_icopyin(mclr,g_nclr,3*mx1);
   gpu_fcopyin(rbufl,g_sbufl,idimp*mcll[3*mx1-1]);
   gpu_fcopyin(rbufr,g_sbufr,idimp*mclr[3*mx1-1]);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[1] += time;
/* second part of particle reorder on x and y cell with mx, my tiles */
   dtimer(&dtime,&itime,-1);
   cgpupppord2lb(g_ppart,g_ppbuff,g_sbufl,g_sbufr,g_kpic,g_ncl,g_ihole,
                 g_ncll,g_nclr,idimp,nppmx,mx1,myp1,npbmx,ntmax,nbmax,
                 g_irc);
   dtimer(&dtime,&itime,1);
   time = (float) dtime;
   ttp[0] += time;
   return;
}

/*--------------------------------------------------------------------*/
void cgpptpose(cuda::std::complex<float> g_bsm[], cuda::std::complex<float> g_btm[],
               cuda::std::complex<float> sm[], cuda::std::complex<float> tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp) {
/* this subroutine sends and receives data between GPUS on different MPI
   nodes to perform a transpose of a matrix distributed in y, to another
   matrix distributed in x.
   one message is sent and received at a time.
   data from GPU is sent asynchronous, overlapping with MPI
   g_bsm/g_btm are complex buffers on GPU to be sent/received
   sm/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   kxyp = kxp*kyp;
/* special case for one processor  */
/* better to use a kernel function */
   if (nvp==1) {
      gpu_ccopyout(sm,g_bsm,kxyp);
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      gpu_ccopyin(tm,g_btm,kxyp);
      return;
   }
   nn = -1;
   ns = kxyp;
   stp = 1; st = 2;
/* send first group to host from GPU */
   gpu_cascopyout(sm,g_bsm,0,ns,stp);
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* send next group to host from GPU */
         if (nn < nvp) {
            st = stp + 1;
            if (st > 2)
               st -= 2;
            gpu_cascopyout(&sm[kxyp*(nn+1)],g_bsm,ns*(nn+1),ns,st);
         }
/* wait for previous MPI sends and receives to complete */
         if (nn > 0) {
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
            cacsndrec(sm,0,0,0,3);
/* copy received group from host to GPU */
            gpu_cascopyin(&tm[kxyp*(nn-1)],g_btm,ns*(nn-1),ns,3);
         }
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
/* ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&mrid); */
         cacsndrec(&tm[kxyp*nn],id,kxyp,n+1,1);
/* wait for previous group to arrive from GPU */
         gpu_waitstream(stp);
         stp = st;
/* send data */
/* ierr = MPI_Isend(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp,&msid); */
         cacsndrec(&sm[kxyp*nn],id,ld,n+1,2);
      }
   }
/* wait for last sends and receives to complete */
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
   cacsndrec(sm,0,0,0,3);
   gpu_cascopyin(&tm[kxyp*nn],g_btm,ns*nn,ns,3);
/* wait for last group item to arrive */
   gpu_waitstream(3);
   return;
}

/*--------------------------------------------------------------------*/
void cgpptposen(cuda::std::complex<float> g_bsm[], cuda::std::complex<float> g_btm[],
                cuda::std::complex<float> sm[], cuda::std::complex<float> tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim) {
/* this subroutine sends and receives data between GPUS on different MPI
   nodes to perform a transpose of an n component  matrix distributed in
   y, to another an n component matrix distributed in x.
   one message is sent and received at a time.
   data from GPU is sent asynchronous, overlapping with MPI.
   g_bsm/g_btm are complex buffers on GPU to be sent/received
   sm/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int  j, n, nn, ks, kyps, kxyp, id, joff, ld, ns, st, stp;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = ndim*(kyp < kyps ? kyp : kyps);
   kxyp = kxp*ndim*kyp;
/* special case for one processor  */
/* better to use a kernel function */
   if (nvp==1) {
      gpu_ccopyout(sm,g_bsm,kxyp);
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      gpu_ccopyin(tm,g_btm,kxyp);
      return;
   }
   nn = -1;
   ns = kxyp;
   stp = 1; st = 2;
/* send first group to host from GPU */
   gpu_cascopyout(sm,g_bsm,0,ns,stp);
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* send next group to host from GPU */
         if (nn < nvp) {
            if (st > 2)
               st -= 2;
            gpu_cascopyout(&sm[kxyp*(nn+1)],g_bsm,ns*(nn+1),
                              ns,st);
         }
/* wait for previous MPI sends and receives to complete */
         if (nn > 0) {
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
            cacsndrec(sm,0,0,0,3);
/* copy received group from host to GPU */
            gpu_cascopyin(&tm[kxyp*(nn-1)],g_btm,ns*(nn-1),
                             ns,3);
         }
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
/* ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&mrid); */
         cacsndrec(&tm[kxyp*nn],id,kxyp,n+1,1);
/* wait for previous group to arrive from GPU */
         gpu_waitstream(stp);
         stp = st;
/* send data */
/* ierr = MPI_Isend(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp,&msid); */
         cacsndrec(&sm[kxyp*nn],id,ld,n+1,2);
      }
   }
/* wait for last sends and receives to complete */
/* ierr = MPI_Wait(&msid,&istatus); ierr = MPI_Wait(&mrid,&istatus); */
   cacsndrec(sm,0,0,0,3);
   gpu_cascopyin(&tm[kxyp*nn],g_btm,ns*nn,ns,3);
/* wait for last group item to arrive */
   gpu_waitstream(3);
   return;
}
