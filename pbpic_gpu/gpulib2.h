/* C header file for CUDA utility Library */

#include <complex.h>

void gpu_setgbsize(int nblock);

int getmmcc();

void gpu_fallocate(float **g_f, int nsize, int *irc);

void gpu_iallocate(int **g_i, int nsize, int *irc);

void gpu_callocate(std::complex<float> **g_c, int nsize, int *irc);

void gpu_deallocate(void *g_d, int *irc);

void hpl_fallocate(float **h_f, int nsize, int *irc);

void hpl_callocate(std::complex<float> **h_c, int nsize, int *irc);

void hpl_deallocate(void *h_d, int *irc);

void gpu_fcopyin(float f[], float g_f[], int nsize);

void gpu_fcopyout(float f[], float g_f[], int nsize);

void gpu_icopyin(int f[], int g_f[], int nsize);

void gpu_icopyout(int f[], int g_f[], int nsize);

void gpu_ccopyin(std::complex<float> f[], std::complex<float> g_f[], int nsize);

void gpu_ccopyout(std::complex<float> f[], std::complex<float> g_f[], int nsize);

void gpu_initstream(int nstream);

void gpu_delstream(int nstream);

void gpu_waitstream(int nstream);

void gpu_cascopyin(std::complex<float> *f, std::complex<float> *g_f, int noff, 
                   int nsize, int nstream);

void gpu_cascopyout(std::complex<float> *f, std::complex<float> *g_f, int noff,
                    int nsize, int nstream);

void gpu_zfmem(float g_f[], int nsize);

void gpu_zcmem(std::complex<float> g_f[], int nsize);

void gpu_set_cache_size(int nscache);

void emptykernel();

void init_cu(int dev, int *irc);

void end_cu();

