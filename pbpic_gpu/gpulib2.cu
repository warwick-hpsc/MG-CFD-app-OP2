/* CUDA utility Library */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
/*#include <complex.h>*/
#include <cuda/std/complex>
#include "cuda.h"

int nblock_size = 64;
int ngrid_size = 1;
int maxgsx = 65535;
int mmcc = 0;
static int devid;

static cudaError_t crc;

#define MAXSTREAMS             4
static cudaStream_t streams[MAXSTREAMS] = {NULL,NULL,NULL,NULL};

__global__ void emptyKernel() {}

/*--------------------------------------------------------------------*/
void gpu_setgbsize(int nblock) {
/* set blocksize */
   nblock_size = nblock;
   return;
}

/*--------------------------------------------------------------------*/
int getmmcc() {
/* get major and minor computer capability */
   return mmcc;
}

/*--------------------------------------------------------------------*/
void gpu_fallocate(float **g_f, int nsize, int *irc) {
/* allocate global float memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMalloc float Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_f = (float *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
void gpu_iallocate(int **g_i, int nsize, int *irc) {
/* allocate global integer memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(int)*nsize);
   if (crc) {
      printf("cudaMalloc int Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_i = (int *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
void gpu_callocate(cuda::std::complex<float> **g_c, int nsize, int *irc) {
/* allocate global float memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(cuda::std::complex<float>)*nsize);
   if (crc) {
      printf("cudaMalloc cuda::std::complex<float> Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_c = (cuda::std::complex<float> *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
void gpu_deallocate(void *g_d, int *irc) {
/* deallocate global memory on GPU */
   crc = cudaFree(g_d);
   if (crc) {
      printf("cudaFree Error=%d:%s\n",crc,cudaGetErrorString(crc));
      *irc = 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void hpl_fallocate(float **h_f, int nsize, int *irc) {
/* allocate page-locked float memory on host, return pointer to C */
   void *hptr = NULL;
   crc = cudaMallocHost(&hptr,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMallocHost float Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *h_f = (float *)hptr;
   return;
}

/*--------------------------------------------------------------------*/
void hpl_callocate(cuda::std::complex<float> **h_c, int nsize, int *irc) {
/* allocate page-locked float memory on host, return pointer to C */
   void *hptr = NULL;
   crc = cudaMallocHost(&hptr,sizeof(cuda::std::complex<float>)*nsize);
   if (crc) {
      printf("cudaMallocHost cuda::std::complex<float> Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *h_c = (cuda::std::complex<float> *)hptr;
   return;
}

/*--------------------------------------------------------------------*/
void hpl_deallocate(void *h_d, int *irc) {
/* deallocate page-locked on host */
   crc = cudaFreeHost(h_d);
   if (crc) {
      printf("cudaFreeHost Error=%d:%s\n",crc,cudaGetErrorString(crc));
      *irc = 1;
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_fcopyin(float *f, float *g_f, int nsize) {
/* copy float array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(float)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice float Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_fcopyout(float *f, float *g_f, int nsize) {
/* copy float array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(float)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost float Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_icopyin(int *f, int *g_f, int nsize) {
/* copy int array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(int)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice int Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_icopyout(int *f, int *g_f, int nsize) {
/* copy int array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(int)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost int Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_ccopyin(cuda::std::complex<float> *f, cuda::std::complex<float> *g_f, int nsize) {
/* copy float array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(cuda::std::complex<float>)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice cuda::std::complex<float> Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_ccopyout(cuda::std::complex<float> *f, cuda::std::complex<float> *g_f, int nsize) {
/* copy cuda::std::complex<float> array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(cuda::std::complex<float>)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost cuda::std::complex<float> Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_initstream(int nstream) {
/* Create Stream for requested identifier nstream       */
/* nstream should be between 1 and MAXSTREAMS inclusive */
   if ((nstream < 1) || (nstream > MAXSTREAMS)) {
      printf("gpu_initstream: nstream out of bounds = %d\n",nstream);
      exit(1);
   }
   if (streams[nstream-1] != NULL) {
      printf("gpu_initstream: nstream already used = %d\n",nstream);
      exit(1);
   }
   crc = cudaStreamCreate(&streams[nstream-1]);
   if (crc) {
      printf("cudaStreamCreate Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_delstream(int nstream) {
/* Destroy Stream for requested identifier nstream      */
/* nstream should be between 1 and MAXSTREAMS inclusive */
   if ((nstream < 1) || (nstream > MAXSTREAMS)) {
      printf("gpu_delstream: nstream out of bounds = %d\n",nstream);
   }
   if (streams[nstream-1] == NULL) {
      printf("gpu_delstream: nstream not allocated = %d\n",nstream);
   }
   crc = cudaStreamDestroy(streams[nstream-1]);
   if (crc) {
      printf("cudaStreamDestroy Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_waitstream(int nstream) {
/* Synchronize Stream for requested identifier nstream  */
/* nstream should be between 0 and MAXSTREAMS inclusive */
   cudaStream_t stream = NULL;
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_waitstream: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaStreamSynchronize(stream);
   if (crc) {
      printf("cudaStreamSynchronize Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_cascopyin(cuda::std::complex<float> *f, cuda::std::complex<float> *g_f, int noff, 
                              int nsize, int nstream) {
/* copy cuda::std::complex<float> array segment from host memory to global GPU memory */
/* asynchronous copy */
   cuda::std::complex<float> *cptr;
   cudaStream_t stream = NULL;
   cptr = &g_f[noff];
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_cascopyin: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaMemcpyAsync((void *)cptr,f,sizeof(cuda::std::complex<float>)*nsize,
                         cudaMemcpyHostToDevice,stream);
   if (crc) {
      printf("Async cudaMemcpyHostToDevice cuda::std::complex<float> Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_cascopyout(cuda::std::complex<float> *f, cuda::std::complex<float> *g_f, int noff,
                               int nsize, int nstream) {
/* copy cuda::std::complex<float> array segment from global GPU memory to host memory */
/* asynchronous copy */
   cuda::std::complex<float> *cptr;
   cudaStream_t stream = NULL;
   cptr = &g_f[noff];
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_cascopyout: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaMemcpyAsync(f,(void *)cptr,sizeof(cuda::std::complex<float>)*nsize,
                         cudaMemcpyDeviceToHost,stream);
   if (crc) {
      printf("Async cudaMemcpyDeviceToHost cuda::std::complex<float> Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_zfmem(float *g_f, int nsize) {
/* initialize float array in global GPU memory to zero */
   crc = cudaMemset((void *)g_f,0,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMemset Error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_zcmem(cuda::std::complex<float> *g_f, int nsize) {
/* initialize cuda::std::complex<float> array in global GPU memory to zero */
   crc = cudaMemset((void *)g_f,0,sizeof(cuda::std::complex<float>)*nsize);
   if (crc) {
      printf("cudaMemset Error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void gpu_set_cache_size(int nscache) {
/* request preferred cache size, requires CUDA 3.2 or higher */
/* nscache = (0,1,2) = (no,small,big) cache size */
   cudaFuncCache cpref;
   if ((nscache < 0) || (nscache > 2))
      return;
   if (nscache==0)
      cpref = cudaFuncCachePreferNone;
   else if (nscache==1)
      cpref = cudaFuncCachePreferShared;
   else if (nscache==2)
      cpref = cudaFuncCachePreferL1;
   crc = cudaThreadSetCacheConfig(cpref);
/* crc = cudaDeviceSetCacheConfig(cpref); */
   if (crc) {
      printf("cudaThreadSetCacheConfig error=%d:%s\n",crc,
             cudaGetErrorString(crc));
   }
   return;
}

/*--------------------------------------------------------------------*/
void emptykernel() {
   int ngx, ngy;
   ngx  = nblock_size < 32768 ? nblock_size : 32768;
   ngy = (ngrid_size - 1)/ngx + 1;
   dim3 dimBlock(nblock_size,1);
   dim3 dimGrid(ngx,ngy);
   crc = cudaGetLastError();
   emptyKernel<<<dimGrid,dimBlock>>>();
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("emptyKernel error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void init_cu(int dev, int *irc, int proc, FILE *fp) {
/* initialize CUDA with device dev or selects best GPU available       */
/* searches throughs devices, selects the device with the most compute */
/* units, and saves the device id devid                                */
/* if dev is a valid device, it is used, otherwise the GPU with the    */
/* most multi-processors is selected                                   */
/* error code is modified only if there is an error */
   int maxcpus = 0, jm = -1;
   int j, ndevs, maxunits;
   unsigned long msize;
   double z;
   struct cudaDeviceProp prop;
/* returns number of device */
   crc = cudaGetDeviceCount(&ndevs);
   if (crc) {
      printf("cudaGetDeviceCount Error=%i:%s\n",crc,
             cudaGetErrorString(crc));
      *irc = 1;
      return;
   }
/* get information about devices */
   for (j = 0; j < ndevs; j++) {
      crc = cudaGetDeviceProperties(&prop,j);
      if (crc) {
         printf("cudaGetDeviceProperties Error=%i:%s\n",crc,
                cudaGetErrorString(crc));
         prop.name[0] = 0;
      }
      maxunits = prop.multiProcessorCount;
      if (dev <= 0) {
         fprintf(fp,"j=%i:CUDA_DEVICE_NAME=%s,CUDA_MULTIPROCESSOR_COUNT=%i\n",
                j,prop.name,maxunits);
         msize = prop.totalGlobalMem;
         z = ((double) msize)/1073741824.0;
         mmcc = 10*prop.major + prop.minor;
         fprintf(fp,"    CUDA_GLOBAL_MEM_SIZE=%lu(%f GB),Capability=%d\n",
                msize,(float) z,mmcc);
         if (maxunits > maxcpus) {
            maxcpus = maxunits;
            jm = j;
         }
      }
   }
   devid = jm;
   if (dev >= 0)
      devid = dev % ndevs;
   fprintf(fp, "proc %i using device j=%i\n",proc, devid);
/* get properties for this device */
   crc = cudaGetDeviceProperties(&prop,devid);
   maxgsx = prop.maxGridSize[0];
   mmcc = 10*prop.major + prop.minor;
/* set device */
   crc = cudaSetDevice(devid);
   if (crc) {
      printf("cudaSetDevice Error=%i:%s\n",crc,
             cudaGetErrorString(crc));
      *irc = 1;
      return;
   }
/* run empty kernel */
   emptykernel();
   return;
}

void end_cu() {
/* terminate CUDA */
   crc = cudaThreadExit();
   if (crc) {
      printf("cudaThreadExit Error=%d:%s\n",crc,cudaGetErrorString(crc));
   }
   return;
}
