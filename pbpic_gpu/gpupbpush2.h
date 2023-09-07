/* C header file for gpupbpush2.cu */

void cgpuppgbppush23l(float ppart[], float fxy[], float bxy[],
                      int kpic[], int noff, int nyp, float qbm,
                      float dt, float dtc, float *ek, int idimp,
                      int nppmx, int nx, int ny, int mx, int my, 
                      int nxv, int nypmx, int mx1, int mxyp1, int ipbc);

void cgpuppgrbppush23l(float ppart[], float fxy[], float bxy[],
                       int kpic[], int noff, int nyp, float qbm,
                       float dt, float dtc, float ci, float *ek,
                       int idimp, int nppmx, int nx, int ny, int mx,
                       int my, int nxv, int nypmx, int mx1, int mxyp1,
                       int ipbc);

void cgpu2ppgppost2l(float ppart[], float q[], int kpic[], int noff,
                     float qm, int idimp, int nppmx, int mx, int my,
                     int nxv, int nypmx, int mx1, int mxyp1);

void cgpu2ppjppost2l(float ppart[], float cu[], int kpic[], int noff,
                     float qm, float dt, int nppmx, int idimp, int nx,
                     int ny, int mx, int my, int nxv, int nypmx,
                     int mx1, int mxyp1, int ipbc);

void cgpu2pprjppost2l(float ppart[], float cu[], int kpic[], int noff,
                      float qm, float dt, float ci, int nppmx,
                      int idimp, int nx, int ny, int mx, int my,
                      int nxv, int nypmx, int mx1, int mxyp1, int ipbc);

void cgpuppcaguard2xl(cuda::std::complex<float> qc[], float scs[], float q[],
                      int nyp, int nx, int nxe, int nypmx, int nxvh,
                      int kypd);

void cgpuppcaguard2yl(cuda::std::complex<float> fc[], float scr[], int nx, int nxvh,
                      int kypd);

void cgpuppcacguard2xl(cuda::std::complex<float> cuc[], float scs[], float cu[],
                       int nyp, int nx, int nxe, int nypmx, int nxvh,
                       int kypd);

void cgpuppcacguard2yl(cuda::std::complex<float> fvc[], float scr[], int nx,
                       int nxvh, int kypd);

void cgpuppcbguard2xl(cuda::std::complex<float> fxyc[], float scs[], float fxy[],
                      int nyp, int nx, int nxe, int nypmx, int nxvh,
                      int kypd);

void cgpuppcbguard2yl(float fxy[], float scr[], int nyp, int nx,
                      int nxe, int nxvh, int nypmx);

void cgpupppord2la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int kpic[], int ncl[], int ihole[], 
                   int ncll[], int nclr[], int noff, int nyp, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int mx1,
                   int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void cgpupppord2lb(float ppart[], float ppbuff[], float rbufl[],
                   float rbufr[], int kpic[], int ncl[], int ihole[],
                   int mcll[], int mclr[], int idimp, int nppmx, int mx1,
                   int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void cgpuppemfield2t(cuda::std::complex<float> fxyt[], cuda::std::complex<float> exyt[],
                     cuda::std::complex<float> ffct[], int isign, int nx, int ny,
                     int kstrt, int nyv, int kxp1, int nyhd);

void cgpuwppfft2rcsx(cuda::std::complex<float> f[], cuda::std::complex<float> bsm[], int isign,
                     int mixup[], cuda::std::complex<float> sct[], int indx,
                     int indy, int kstrt, int nvp, int kxp1, int kyp,
                     int nxhd, int kypd, int nxhyd, int nxyhd);

void cgpuwppfft2rcsy(cuda::std::complex<float> g[], cuda::std::complex<float> brm[], int isign,
                     int mixup[], cuda::std::complex<float> sct[], int indx,
                     int indy, int kstrt, int nvp, int kxp1, int kyp,
                     int nyd, int nxhyd, int nxyhd);

void cgpuwppfft2rcsxn(cuda::std::complex<float> fn[], cuda::std::complex<float> bsm[], int isign,
                      int mixup[], cuda::std::complex<float> sct[], int indx,
                      int indy, int ndim, int kstrt, int nvp, int kxp1,
                      int kyp, int nxhd, int kypd, int nxhyd, int nxyhd);

void cgpuwppfft2rcsyn(cuda::std::complex<float> gn[], cuda::std::complex<float> brm[], int isign,
                      int mixup[], cuda::std::complex<float> sct[], int indx,
                      int indy, int ndim, int kstrt, int nvp, int kxp1,
                      int kyp, int nyd, int nxhyd, int nxyhd);

void cgpuppltpose(cuda::std::complex<float> f[], cuda::std::complex<float> g[], int nx, int ny,
                  int kxp, int kyp, int kstrt, int nxv, int nyv);

void cgpuppltposen(cuda::std::complex<float> fn[], cuda::std::complex<float> gn[], int nx,
                   int ny, int kxp, int kyp, int kstrt, int ndim,
                   int nxv, int nyv);

void cgpusum2(float a[], float *sa, int nx);

