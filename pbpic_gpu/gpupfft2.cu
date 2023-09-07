/* CUDA Parallel FFT Library */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"
#include <cuda/std/complex>
#include <cufft.h>

extern int nblock_size;
extern int maxgsx;

static cudaError_t crc;
static cufftResult cfrc = CUFFT_SUCCESS;
static cufftHandle planrx = 0, planxr = 0, planrxn = 0, planxrn = 0;
static cufftHandle plany = 0, planyn = 0;


__global__ void gpuppmtposes(cuda::std::complex<float> f[], cuda::std::complex<float> sm[], int nx, int kxp,
                             int kyps, int kstrt, int nvp, int kxyp,
                             int nxv, int kypd);

__global__ void gpuppmtposer(cuda::std::complex<float> g[], cuda::std::complex<float> tm[], int ny, int kyp,
                             int kxps, int kstrt, int nvp, int kxyp,
                             int nyv, int kxpd);
                             
__global__ void gpuppmtposesn(cuda::std::complex<float> fn[], cuda::std::complex<float> sm[], int nx, int kxp,
                              int kyps, int kstrt, int nvp, int ndim,
                              int kxyp, int nxv, int kypd);

__global__ void gpuppmtposern(cuda::std::complex<float> gn[], cuda::std::complex<float> tm[], int ny, int kyp,
                              int kxps, int kstrt, int nvp, int ndim,
                              int kxyp, int nyv, int kxpd);


/*--------------------------------------------------------------------*/
__global__ void gpuppsmtposes(cuda::std::complex<float> f[], cuda::std::complex<float> sm[], float ani,
                              int nx, int kxp, int kyps, int kstrt,
                              int nvp, int kxyp, int nxv, int kypd) {
/* extract data to send and normalize */
/* local data */
   int ks, j, k, n, nn, id, joff, ld;
   cuda::std::complex<float> a;
   ks = kstrt - 1;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.y;
   if (n < nvp) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*ks;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not send local data */
      if (id != ks) {
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kxp < ld ? kxp : ld;
/*       for (k = 0; k < kyps; k++) { */
         k = blockIdx.x;
         if (k < kyps) {
/*          for (j = 0; j < ld; j++) { */
            j = threadIdx.x;
            while (j < ld) {
               a = f[j+joff+nxv*k];
               a.real(ani*a.real());
               a.imag(ani*a.imag());
               sm[j+ld*k+kxyp*n] = a;
               j += blockDim.x;
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppsmtposesn(cuda::std::complex<float> fn[], cuda::std::complex<float> sm[], float ani,
                               int nx, int kxp, int kyps, int kstrt,
                               int nvp, int ndim, int kxyp, int nxv,
                               int kypd) {
/* extract vector data to send and normalize */
/* local data */
   int ks, i, j, k, n, nn, id, joff, ld, nnxv, nkxyp;
   cuda::std::complex<float> a;
   ks = kstrt - 1;
   nnxv = ndim*nxv;
   nkxyp = ndim*kxyp;
/* for (n = 0; n < nvp; n++) { */
   n = blockIdx.y;
   if (n < nvp) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* find which node sends to itself */
      nn = 2*ks;
      if (nn >= nvp)
         nn -= nvp;
/* adjust counter */
      if (n > nn)
         n -= 1;
/* do not send local data */
      if (id != ks) {
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kxp < ld ? kxp : ld;
/*       for (k = 0; k < kyps; k++) { */
         k = blockIdx.x;
         if (k < kyps) {
/*          for (j = 0; j < ld; j++) { */
            j = threadIdx.x;
            while (j < ld) {
               for (i = 0; i < ndim; i++) {
                  a = fn[j+joff+nxv*i+nnxv*k];
                  a.real(ani*a.real());
                  a.imag(ani*a.imag());
                  sm[j+ld*(i+ndim*k)+nkxyp*n] = a;
               }
               j += blockDim.x;
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppsltpose(cuda::std::complex<float> f[], cuda::std::complex<float> g[], float ani, int nx,
                             int ny, int kxp, int kyp, int kstrt,
                             int nxv, int nyv) {
/* transpose local data with scaling */
/* local data */
   int mxv, j, k, ks, kxps, kyps, joff, koff, js, jj, kk;
   cuda::std::complex<float> a;
/* The size of the shared memory array is as follows: */
/* cuda::std::complex<float> s2[(mx + 1)*mx];                            */
   extern __shared__ cuda::std::complex<float> s2[];
   mxv = blockDim.x + 1;
   ks = kstrt - 1;
   joff = kxp*ks;
   koff = kyp*ks;
   kxps = nx - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - koff;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   js = threadIdx.x;
   ks = threadIdx.y;
   jj = blockDim.x*blockIdx.x;
   kk = blockDim.y*blockIdx.y;
   j = js + jj;
   k = ks + kk;
   if ((j < kxps) && (k < kyps)) {
      s2[js+mxv*ks] = f[j+joff+nxv*k];
   }
/* synchronize threads */
   __syncthreads();
   j = ks + jj;
   k = js + kk;
   if ((j < kxps) && (k < kyps)) {
      a = s2[ks+mxv*js];
      a.real(ani*a.real());
      a.imag(ani*a.imag());
      g[k+koff+nyv*j] = a;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpuppsltposen(cuda::std::complex<float> fn[], cuda::std::complex<float> gn[], float ani,
                              int nx, int ny, int kxp, int kyp,
                              int kstrt, int ndim, int nxv, int nyv) {
/* transpose local vector data with scaling */
/* local data */
   int mxv, i, j, k, ks, kxps, kyps, joff, koff, js, jj, kk;
   int nnxv, nnyv;
   cuda::std::complex<float> a;
/* The size of the shared memory array is as follows: */
/* cuda::std::complex<float> s2n[ndim*(mx + 1)*mx];                      */
   extern __shared__ cuda::std::complex<float> s2n[];
   mxv = blockDim.x + 1;
   ks = kstrt - 1;
   nnxv = ndim*nxv;
   nnyv = ndim*nyv;
   joff = kxp*ks;
   koff = kyp*ks;
   kxps = nx - joff;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - koff;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   js = threadIdx.x;
   ks = threadIdx.y;
   jj = blockDim.x*blockIdx.x;
   kk = blockDim.y*blockIdx.y;
   j = js + jj;
   k = ks + kk;
   if ((j < kxps) && (k < kyps)) {
      for (i = 0; i < ndim; i++) {
         s2n[js+mxv*(i+ndim*ks)] = fn[j+joff+nxv*i+nnxv*k];
      }
   }
/* synchronize threads */
   __syncthreads();
   j = ks + jj;
   k = js + kk;
   if ((j < kxps) && (k < kyps)) {
      for (i = 0; i < ndim; i++) {
         a = s2n[ks+mxv*(i+ndim*js)];
         a.real(ani*a.real());
         a.imag(ani*a.imag());
         gn[k+koff+nyv*i+nnyv*j] = a;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/	
void gpupfft2rrcuinit(int nx, int kypp, int ndim) {
   if (kypp <= 0)
      return;
   cfrc = cufftPlan1d(&planrx,nx,CUFFT_R2C,kypp);
   if (cfrc) {
      printf("cufftPlan1d planrx error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planxr,nx,CUFFT_C2R,kypp);
   if (cfrc) {
      printf("cufftPlan1d planxr error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planrxn,nx,CUFFT_R2C,ndim*kypp);
   if (cfrc) {
      printf("cufftPlan1d planrxn error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planxrn,nx,CUFFT_C2R,ndim*kypp);
   if (cfrc) {
      printf("cufftPlan1d planxrn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/	
void gpupfft2cuinit(int kxpp, int ny, int ndim) {
   if (kxpp <= 0)
      return;
   cfrc = cufftPlan1d(&plany,ny,CUFFT_C2C,kxpp);
   if (cfrc) {
      printf("cufftPlan1d plany error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planyn,ny,CUFFT_C2C,ndim*kxpp);
   if (cfrc) {
      printf("cufftPlan1d planyn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2rrcudel() {
   if (planrx != 0)
      cfrc = cufftDestroy(planrx);
   if (cfrc) {
      printf("cufftPlan1d planrx error=%d\n",cfrc);
      exit(1);
   }
   if (planxr != 0)
      cfrc = cufftDestroy(planxr);
   if (cfrc) {
      printf("cufftPlan1d planxr error=%d\n",cfrc);
      exit(1);
   }
   if (planrxn != 0)
      cfrc = cufftDestroy(planrxn);
   if (cfrc) {
      printf("cufftPlan1d planrxn error=%d\n",cfrc);
      exit(1);
   }
   if (planxr != 0)
      cfrc = cufftDestroy(planxrn);
   if (cfrc) {
      printf("cufftPlan1d planxrn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2cudel() {
   if (plany != 0)
      cfrc = cufftDestroy(plany);
   if (cfrc) {
      printf("cufftPlan1d plany error=%d\n",cfrc);
      exit(1);
   }
   if (planyn != 0)
      cfrc = cufftDestroy(planyn);
   if (cfrc) {
      printf("cufftPlan1d planyn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2rrcux(cuda::std::complex<float> *f, cuda::std::complex<float> *bsm, int isign, 
                              int indx, int indy, int kstrt, int nvp,
                              int kxp1, int kyp, int nxh1d, int kypd) {
/* wrapper function for real to complex fft in x,                */
/* without packed data                                           */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* nxh1d must be = nx/2+1                                        */
/* local data */
   int nx, nxh1, ny, ks, kypp, kxyp, ns;
   int mx = 16;
   float ani;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
   if (kypp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGrids(kypp,nvp);
   dim3 dimGridty((kyp-1)/mx+1,(kxp1-1)/mx+1,nvp);
   ns = (mx+1)*mx*sizeof(cuda::std::complex<float>);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfrc = cufftExecR2C(planrx,(cufftReal *)f,(cufftComplex *)f);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecR2C(-1) planrx error=%d\n",cfrc);
         exit(1);
      }
/* extract data to send and normalize */
      ani = 1.0f/(((float) nx)*((float) ny));
      crc = cudaGetLastError();
      gpuppsmtposes<<<dimGrids,dimBlock>>>(f,bsm,ani,nxh1,kxp1,kypp,
                                           kstrt,nvp,kxyp,nxh1d,kypd);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppsmtposes error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposer<<<dimGridty,dimBlockt,ns>>>(f,bsm,nxh1,kxp1,kypp,
                                               kstrt,nvp,kxyp,nxh1d,
                                               kypd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposer error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      cfrc = cufftExecC2R(planxr,(cufftComplex *)f,(cufftReal *)f);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2R(1) planxr error=%d\n",cfrc);
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2rrcuy(cuda::std::complex<float> *g, cuda::std::complex<float> *brm, int isign,
                              int indx, int indy, int kstrt, int nvp, 
                              int kxp1, int kyp, int nyd) {
/* wrapper function for real to complex fft in y,                */
/* without packed data                                           */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* local data */
   int nx, nxh1, ny, ks, kxpp, kxyp, ns;
   int mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh1 - kxp1*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp1 < kxpp ? kxp1 : kxpp;
   if (kxpp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGrids(kxpp,nvp);
   dim3 dimGridtx((kxp1-1)/mx+1,(kyp-1)/mx+1,nvp);
   ns = (mx+1)*mx*sizeof(cuda::std::complex<float>);
/* inverse fourier transform */
   if (isign < 0) {
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposer<<<dimGridtx,dimBlockt,ns>>>(g,brm,ny,kyp,kxpp,kstrt,
                                               nvp,kxyp,nyd,kxp1);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposer error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      cfrc = cufftExecC2C(plany,(cufftComplex *)g,(cufftComplex *)g,
                          CUFFT_FORWARD);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2C(-1) plany error=%d\n",cfrc);
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfrc = cufftExecC2C(plany,(cufftComplex *)g,(cufftComplex *)g,
                          CUFFT_INVERSE);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecC2C(1) plany error=%d\n",cfrc);
         exit(1);
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposes<<<dimGrids,dimBlock>>>(g,brm,ny,kyp,kxpp,kstrt,
                                          nvp,kxyp,nyd,kxp1);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposes error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2rrcuxn(cuda::std::complex<float> *fn, cuda::std::complex<float> *bsm, int isign, 
                               int indx, int indy, int ndim, int kstrt,
                               int nvp, int kxp1, int kyp, int nxh1d,
                               int kypd) {
/* wrapper function for real to complex fft in x,                */
/* without packed data                                           */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* ndim = vector dimension                                       */
/* nxh1d must be = nx/2+1                                        */
/* local data */
   int nx, nxh1, ny, ks, kypp, kxyp, ns;
   int mx = 16;
   float ani;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kypp = ny - kyp*ks;
   kypp = 0 > kypp ? 0 : kypp;
   kypp = kyp < kypp ? kyp : kypp;
   if (kypp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGrids(kypp,nvp);
   dim3 dimGridty((kyp-1)/mx+1,(kxp1-1)/mx+1,nvp);
   ns = ndim*(mx+1)*mx*sizeof(cuda::std::complex<float>);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfrc = cufftExecR2C(planrxn,(cufftReal *)fn,(cufftComplex *)fn);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecR2C(-1) planrxn error=%d\n",cfrc);
         exit(1);
      }
/* extract data to send and normalize */
      ani = 1.0f/(((float) nx)*((float) ny));
      crc = cudaGetLastError();
      gpuppsmtposesn<<<dimGrids,dimBlock>>>(fn,bsm,ani,nxh1,kxp1,kypp,
                                            kstrt,nvp,ndim,kxyp,nxh1d,
                                            kypd);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppsmtposesn error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposern<<<dimGridty,dimBlockt,ns>>>(fn,bsm,nxh1,kxp1,kypp,
                                                kstrt,nvp,ndim,kxyp,
                                                nxh1d,kypd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposern error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      cfrc = cufftExecC2R(planxrn,(cufftComplex *)fn,(cufftReal *)fn);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2R(1) planxrn error=%d\n",cfrc);
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpupfft2rrcuyn(cuda::std::complex<float> *gn, cuda::std::complex<float> *brm, int isign,
                               int indx, int indy, int ndim, int kstrt,
                               int nvp, int kxp1, int kyp, int nyd) {
/* wrapper function for real to complex fft in y,                */
/* without packed data                                           */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* ndim = vector dimension                                       */
/* local data */
   int nx, nxh1, ny, ks, kxpp, kxyp, ns;
   int mx = 16;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   ks = kstrt - 1;
   kxpp = nxh1 - kxp1*ks;
   kxpp = 0 > kxpp ? 0 : kxpp;
   kxpp = kxp1 < kxpp ? kxp1 : kxpp;
   if (kxpp <= 0)
      return;
   kxyp = kxp1*kyp;
   dim3 dimGrids(kxpp,nvp);
   dim3 dimGridtx((kxp1-1)/mx+1,(kyp-1)/mx+1,nvp);
   ns = ndim*(mx+1)*mx*sizeof(cuda::std::complex<float>);
/* inverse fourier transform */
   if (isign < 0) {
/* transpose data received */
      crc = cudaGetLastError();
      gpuppmtposern<<<dimGridtx,dimBlockt,ns>>>(gn,brm,ny,kyp,kxpp,kstrt,
                                                nvp,ndim,kxyp,nyd,kxp1);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposern error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      cfrc = cufftExecC2C(planyn,(cufftComplex *)gn,(cufftComplex *)gn,
                          CUFFT_FORWARD);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2C(-1) planyn error=%d\n",cfrc);
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfrc = cufftExecC2C(planyn,(cufftComplex *)gn,(cufftComplex *)gn,
                          CUFFT_INVERSE);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecC2C(1) planyn error=%d\n",cfrc);
         exit(1);
      }
/* extract data to send */
      crc = cudaGetLastError();
      gpuppmtposesn<<<dimGrids,dimBlock>>>(gn,brm,ny,kyp,kxpp,kstrt,
                                           nvp,ndim,kxyp,nyd,kxp1);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuppmtposesn error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgpuppsltpose(cuda::std::complex<float> *f, cuda::std::complex<float> *g, float ani, int nx,
                              int ny, int kxp, int kyp, int kstrt,
                              int nxv, int nyv) {
/* local complex transpose with scaling */
/* input = f, output = g                */
/* local data */
   int ns;
   static int mx = 16;
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   dim3 dimGridtx((kxp-1)/mx+1,(kyp-1)/mx+1);
   ns = (mx+1)*mx*sizeof(cuda::std::complex<float>);
/* local transpose f to g */
   crc = cudaGetLastError();
   gpuppsltpose<<<dimGridtx,dimBlockt,ns>>>(f,g,ani,nx,ny,kxp,kyp,kstrt,
                                            nxv,nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppsltpose error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cgpuppsltposen(cuda::std::complex<float> *fn, cuda::std::complex<float> *gn, float ani,
                               int nx, int ny, int kxp, int kyp,
                               int kstrt, int ndim, int nxv, int nyv) {
/* local complex vector transpose with scaling */
/* input = fn, output = gn                     */
/* local data */
   int ns;
   static int mx = 16;
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   dim3 dimGridtx((kxp-1)/mx+1,(kyp-1)/mx+1);
   ns = ndim*(mx+1)*mx*sizeof(cuda::std::complex<float>);
/* local transpose f to g */
   crc = cudaGetLastError();
   gpuppsltposen<<<dimGridtx,dimBlockt,ns>>>(fn,gn,ani,nx,ny,kxp,kyp,
                                             kstrt,ndim,nxv,nyv);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gpuppsltposen error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}
