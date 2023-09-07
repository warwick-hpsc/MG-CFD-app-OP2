/*--------------------------------------------------------------------*/
/* Skeleton 2-1/2D Electromagnetic MPI PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <chrono>
#include <sys/time.h>
#include <stdbool.h>
#include "pbpush2.h"
#include "pplib2.h"
#include "cpx_utils.h"
#include "mpi.h"
#include "../src_op/coupler_config.h"
#include "../src/structures.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

void output_pic(int ntime, int nvp, float wt, float tdpost, 
				float tdjpost, float tguard, float tfield, float tfft[], 
				float tpush, float tmov, float tsort, float time,
                int kstrt, FILE *fp, double comp_sec, double tot_sec);

int main(int argc, char *argv[]) {
   //replace with coupler_config/inputs 
   MPI_Comm pic_comm = MPI_COMM_WORLD;
   //open a file to output to.
   //char filename[2];
   char default_name[24] = "PIC_output_instance_";
   //sprintf(filename, "%d", instance_number);
   //strcat(default_name, filename);
   FILE *fp = fopen(default_name, "w");
/* indx/indy = exponent which determines grid points in x/y direction: */
/* nx = 2**indx, ny = 2**indy */
   int indx, indy;
/* npx/npy = number of electrons distributed in x/y direction */
   int npx, npy;
   read_inputs(&indx, &indy, &npx, &npy);
/* ndim = number of velocity coordinates = 3 */
   int ndim = 3;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.04, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction */
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   float vtz = 1.0, vz0 = 0.0;
/* ax/ay = smoothed particle size in x/y direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ay = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 5 */
/* ipbc = particle boundary condition: 1 = periodic */
/* sortime = number of time steps between standard electron sorting */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 5, ipbc = 1, sortime = 50, relativity = 1;
/* idps = number of partition boundaries */
   int idps = 2;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy;
   int ny1, ntime, isign, ierr;
   float qbme, affp, dth;
   double np;

/* declare scalars for MPI code */
   int ntpose = 1;
   int nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn;
   int nyp, noff, npp, nps, nbmax, ntmax;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* cue = electron current density with guard cells */
/* fxyze/bxyze = smoothed electric/magnetic field with guard cells */
   float *cue = NULL, *fxyze = NULL, *bxyze = NULL;
/* exyz/bxyz = transverse electric/magnetic field in fourier space */
   std::complex<float> *exyz = NULL, *bxyz = NULL;
/* qt = scalar charge density field array in fourier space */
   std::complex<float> *qt = NULL;
/* cut = vector current density field array in fourier space */
/* fxyt/bxyt = vector electric/magnetic field in fourier space */
   std::complex<float> *cut = NULL, *fxyt = NULL, *bxyt = NULL;
/* ffc = form factor array for poisson solver */
   std::complex<float> *ffc = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* sct = sine/cosine table for FFT */
   std::complex<float> *sct = NULL;
/* ihole = location of hole left in particle arrays */
   int *ihole = NULL;
/* npic = scratch array for reordering particles */
   int *npic = NULL;
   double wtot[7], work[7];
   int info[7];

/* declare arrays for MPI code: */
/* bs/br = complex send/receive buffers for data transpose */
   std::complex<float> *bs = NULL, *br = NULL;
/* sbufl/sbufr = particle buffers sent to nearby processors */
/* rbufl/rbufr = particle buffers received from nearby processors */
   float *sbufl = NULL, *sbufr = NULL, *rbufl = NULL, *rbufr = NULL;
/* edges[0:1] = lower:upper y boundaries of particle partition */
   float *edges = NULL;
/* scr = guard cell buffer received from nearby processors */
   float *scr = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tmov = 0.0;
   float tfft[2] = {0.0,0.0};
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
   np =  (double) npx*(double) npy;
/* nx/ny = number of grid points in x/y direction */
   nx = 1L<<indx; ny = 1L<<indy;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2;
   nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe;
   nxyh = (nx > ny ? nx : ny)/2; nxhy = nxh > ny ? nxh : ny;
   ny1 = ny + 1;
   qbme = qme;
   affp = (double) nx*(double) ny/np;
   dth = 0.0;

/* nvp = number of MPI ranks */
/* initialize for distributed memory parallel processing */
   cppinit2(&idproc,&nvp,argc,argv,pic_comm);
   kstrt = idproc + 1;
/* check if too many processors */
   if (nvp > ny) {
      if (kstrt==1) {
         fprintf(fp, "Too many processors requested: ny, nvp=%d,%d\n",ny,nvp);
      }
	  exit(1);
   }

/* initialize data for MPI code */
   edges = (float *) malloc(idps*sizeof(float));
/* calculate partition variables: edges, nyp, noff, nypmx              */
/* edges[0:1] = lower:upper boundary of particle partition             */
/* nyp = number of primary (complete) gridpoints in particle partition */
/* noff = lowermost global gridpoint in particle partition             */
/* nypmx = maximum size of particle partition, including guard cells   */
/* nypmn = minimum value of nyp                                        */
   cpdicomp2l(edges,&nyp,&noff,&nypmx,&nypmn,ny,kstrt,nvp,idps);
   if (nypmn < 1) {
      if (kstrt==1) {
        fprintf(fp, "combination not supported nvp, ny = %d,%d\n",nvp,ny);
      }
	  exit(1);
   }
/* initialize additional scalars for MPI code */
/* kxp = number of complex grids in each field partition in x direction */
   kxp = (nxh - 1)/nvp + 1;
/* kyp = number of complex grids in each field partition in y direction */
   kyp = (ny - 1)/nvp + 1;
/* npmax = maximum number of electrons in each partition */
   npmax = (np/nvp)*1.25;
/* nbmax = size of buffer for passing particles between processors */
   nbmax = 0.1*npmax;
/* ntmax = size of ihole buffer for particles leaving processor */
   ntmax = 2*nbmax;

/* allocate data for standard code */
   part = (float *) malloc(idimp*npmax*sizeof(float));
   if (sortime > 0)
      part2 = (float *) malloc(idimp*npmax*sizeof(float));
   qe = (float *) malloc(nxe*nypmx*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   exyz = (std::complex<float> *) malloc(ndim*nye*kxp*sizeof(std::complex<float>));
   bxyz = (std::complex<float> *) malloc(ndim*nye*kxp*sizeof(std::complex<float>));
   qt = (std::complex<float> *) malloc(nye*kxp*sizeof(std::complex<float>));
   fxyt = (std::complex<float> *) malloc(ndim*nye*kxp*sizeof(std::complex<float>));
   cut = (std::complex<float> *) malloc(ndim*nye*kxp*sizeof(std::complex<float>));
   bxyt = (std::complex<float> *) malloc(ndim*nye*kxp*sizeof(std::complex<float>));
   ffc = (std::complex<float> *) malloc(nyh*kxp*sizeof(std::complex<float>));
   mixup = (int *) malloc(nxhy*sizeof(int));
   sct = (std::complex<float> *) malloc(nxyh*sizeof(std::complex<float>));
   ihole = (int *) malloc((ntmax+1)*sizeof(int));
   npic = (int *) malloc(nypmx*sizeof(int));

/* allocate data for MPI code */
   bs = (std::complex<float> *) malloc(ndim*kxp*kyp*sizeof(std::complex<float>));
   br = (std::complex<float> *) malloc(ndim*kxp*kyp*sizeof(std::complex<float>));
   sbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   sbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   scr = (float *) malloc(ndim*nxe*sizeof(float));

/* prepare fft tables */
   cwpfft2rinit(mixup,sct,indx,indy,nxhy,nxyh);
/* calculate form factors */
   isign = 0;
   cppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,nyh);
/* initialize electrons */
   nps = 1;
   npp = 0;
   cpdistr2h(part,edges,&npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,nx,ny,
             idimp,npmax,idps,ipbc,&ierr);
/* check for particle initialization error */
   if (ierr != 0) {
      if (kstrt==1) {
		fprintf(fp, "particle initialization error: ierr=%d\n",ierr);
      }
	  exit(1);
   }

/* initialize transverse electromagnetic fields */
   for (j = 0; j < ndim*nye*kxp; j++) {
      exyz[j] = std::complex<double>(0.0, 0.0);
      bxyz[j] = std::complex<double>(0.0, 0.0);
   }

   if (dt > 0.45*ci) {
      if (kstrt==1) {
		fprintf(fp, "Warning: Courant condition may be exceeded!\n");
      }
   }

/* * * * start main iteration loop * * * */
   int number_loops = coupler_cycles * pic_conversion_factor;
   std::chrono::duration<double> comp_seconds;
   std::chrono::duration<double> total_seconds;
   std::chrono::steady_clock::time_point start_comp;
   std::chrono::steady_clock::time_point end_comp;
   std::chrono::steady_clock::time_point start_total;
   std::chrono::steady_clock::time_point end_total;
   start_total = std::chrono::steady_clock::now();
   for(ntime = 0; ntime < number_loops; ntime++){
/* output the current iteration and how many we are doing */
	  if(kstrt==1){
		fprintf(fp, "performing iteration %d of %d\n", ntime+1, number_loops);
	  }
	  start_comp = std::chrono::steady_clock::now();
/* deposit current with standard procedure: updates part, cue, ihole */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx; j++) {
         cue[j] = 0.0;
      }
      if (relativity==1) {
         cppgrjpost2l(part,cue,edges,npp,noff,ihole,qme,dth,ci,nx,ny,
                      idimp,npmax,nxe,nypmx,idps,ntmax,ipbc);
      }
      else {
         cppgjpost2l(part,cue,edges,npp,noff,ihole,qme,dth,nx,ny,idimp,
                     npmax,nxe,nypmx,idps,ntmax,ipbc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;
/* check for ihole overflow error */
      if (ihole[0] < 0) {
         ierr = -ihole[0];
		 fprintf(fp, "ihole overflow error: ntmax,ih=%d,%d\n",ntmax,ierr);
         exit(1);
      }
/* move electrons into appropriate spatial regions: updates part, npp */
      dtimer(&dtime,&itime,-1);
      cppmove2(part,edges,&npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp,
               idimp,npmax,idps,nbmax,ntmax,info,fp,debug);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* check for particle manager error */
      if (info[0] != 0) {
         ierr = info[0];
         if (kstrt==1) {
            fprintf(fp, "current particle manager error: ierr=%d\n",ierr);
         }
		 exit(1);
      }

/* deposit charge with standard procedure: updates qe */
      for(int i = 0; i < 12; i++){
	     dtimer(&dtime,&itime,-1);
         for (j = 0; j < nxe*nypmx; j++) {
            qe[j] = 0.0;
         }
         cppgpost2l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tdpost += time;
	  }

/* add guard cells with standard procedure: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      cppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx);
      cppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
      cppaguard2xl(qe,nyp,nx,nxe,nypmx);
      cppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;
/* transform charge to fourier space with standard procedure: updates qt */
/* modifies qe */
      /*dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2r((std::complex<float> *)qe,qt,bs,br,isign,ntpose,mixup,sct,&ttp,
                indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;*/

/* transform current to fourier space with standard procedure: updates cut */
/* modifies cue */
      /*dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2r3((std::complex<float> *)cue,cut,bs,br,isign,ntpose,mixup,sct,
                 &ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,
                 nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;*/

/* take transverse part of current with standard procedure: updates cut */
/* modifies cut */
      /*dtimer(&dtime,&itime,-1);
      cppcuperp2(cut,nx,ny,kstrt,nye,kxp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;*/

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates exyz, bxyz, wf, wm */
      /*dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         cippbpoisp23(cut,bxyz,ffc,ci,&wm,nx,ny,kstrt,nye,kxp,nyh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         cppmaxwel2(exyz,bxyz,cut,ffc,affp,ci,dt,&wf,&wm,nx,ny,kstrt,
                    nye,kxp,nyh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;*/

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyt, we */
      /*dtimer(&dtime,&itime,-1);
      isign = -1;
      cppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,
               nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;*/

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cppemfield2(fxyt,exyz,ffc,isign,nx,ny,kstrt,nye,kxp,nyh);
/* copy magnetic field with standard procedure: updates bxyt */
      isign = -1;
      cppemfield2(bxyt,bxyz,ffc,isign,nx,ny,kstrt,nye,kxp,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with standard procedure: updates fxyze */
/* modifies fxyt */
      /*dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2r3((std::complex<float> *)fxyze,fxyt,bs,br,isign,ntpose,mixup,
                 sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                 nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;*/

/* transform magnetic field to real space with standard procedure: */
/* updates bxyze, modifies bxyt */
      /*dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2r3((std::complex<float> *)bxyze,bxyt,bs,br,isign,ntpose,mixup,
                 sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                 nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;*/

/* copy guard cells with standard procedure: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      cppncguard2l(fxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(fxyze,nyp,nx,ndim,nxe,nypmx);
      cppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles: updates part, wke, and ihole */
      for(int i = 0; i < 30; i++){
	     dtimer(&dtime,&itime,-1);
         wke = 0.0;
         if (relativity==1) {
            cppgrbpush23l(part,fxyze,bxyze,edges,npp,noff,ihole,qbme,dt,
                          dth,ci,&wke,nx,ny,idimp,npmax,nxe,nypmx,idps,
                          ntmax,ipbc);
         }
         else {
            cppgbpush23l(part,fxyze,bxyze,edges,npp,noff,ihole,qbme,dt,dth,
                         &wke,nx,ny,idimp,npmax,nxe,nypmx,idps,ntmax,ipbc);
         }
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tpush += time;
	  }
/* check for ihole overflow error */
      if (ihole[0] < 0) {
         ierr = -ihole[0];
         fprintf(fp, "ihole overflow error: ntmax,ih=%d,%d\n",ntmax,ierr);
         exit(1);
      }

/* move electrons into appropriate spatial regions: updates part, npp */
      dtimer(&dtime,&itime,-1);
      cppmove2(part,edges,&npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp,
               idimp,npmax,idps,nbmax,ntmax,info,fp,debug);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* check for particle manager error */
      if (info[0] != 0) {
         ierr = info[0];
         if (kstrt==1) {
            fprintf(fp, "push particle manager error: ierr=%d\n",ierr);
         }
		 exit(1);
      }

/* sort particles for standard code: updates part */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            cppdsortp2yl(part,part2,npic,npp,noff,nyp,idimp,npmax,nypmx);
/* exchange pointers */
            tpart = part;
            part = part2;
            part2 = tpart;
            dtimer(&dtime,&itime,1);
            time = (float) dtime;
            tsort += time;
         }
      }

/* energy diagnostic */
      wt = we + wf + wm;
      wtot[0] = wt;
      wtot[1] = wke;
      wtot[2] = 0.0;
      wtot[3] = wke + wt;
      wtot[4] = we;
      wtot[5] = wf;
      wtot[6] = wm;
      cppdsum(wtot,work,7);
      wke = wtot[1];
      we = wtot[4];
      wf = wtot[5];
      wm = wtot[6];
   
      end_comp = std::chrono::steady_clock::now();
      comp_seconds += (end_comp - start_comp);
   }
   end_total = std::chrono::steady_clock::now();
   total_seconds = end_total - start_total;
   output_pic(ntime, nvp, wt, tdpost, tdjpost, tguard, 
			  tfield, tfft, tpush, tmov, tsort, time, 
			  kstrt, fp, comp_seconds.count(), 
			  total_seconds.count());
}
/* * * * end main iteration loop * * * */
void output_pic(int ntime, int nvp, float wt, float tdpost, 
				float tdjpost, float tguard, float tfield, float tfft[], 
				float tpush, float tmov, float tsort, float time,
				int kstrt, FILE *fp, double comp_sec, double tot_sec){
   if (kstrt==1) {
	 fprintf(fp, "\n\n\n-----------------------------\n");
     fprintf(fp, "Number of cycles = %i \n", ntime);
     fprintf(fp, "MPI nodes = %i\n", nvp); 

	 tdpost += tdjpost;
	 tfield += tguard + tfft[0];
	 tsort += tmov;
	 time = tdpost + tpush + tsort;
	 wt = time + tfield;
	 fprintf(fp, "total time = %f\n",wt);
	 fprintf(fp, "total time in loop = %f\n", tot_sec);
	 fprintf(fp, "total comp time in loop = %f\n", comp_sec);
   }
}
