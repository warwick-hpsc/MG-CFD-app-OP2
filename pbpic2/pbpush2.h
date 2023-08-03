/* header file for pbpush2.c */

double ranorm();

void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps);

void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr);

void cppgbpush23l(float part[], float fxy[], float bxy[], float edges[],
                  int npp, int noff, int ihole[], float qbm, float dt, 
                  float dtc, float *ek, int nx, int ny, int idimp,
                  int npmax, int nxv, int nypmx, int idps, int ntmax,
                  int ipbc);

void cdppgbpush23l(double part[], double fxy[], double bxy[], 
                   double edges[], int npp, int noff, int ihole[],
                   double qbm, double dt, double dtc, double *ek,
                   int nx, int ny, int idimp, int npmax, int nxv,
                   int nypmx, int idps, int ntmax, int ipbc);

void cppgrbpush23l(float part[], float fxy[], float bxy[],
                   float edges[], int npp, int noff, int ihole[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int nx, int ny, int idimp, int npmax, int nxv,
                   int nypmx, int idps, int ntmax, int ipbc);

void cdppgrbpush23l(double part[], double fxy[], double bxy[],
                    double edges[], int npp, int noff, int ihole[],
                    double qbm, double dt, double dtc, double ci, 
                    double *ek, int nx, int ny, int idimp, int npmax,
                    int nxv, int nypmx, int idps, int ntmax, int ipbc);

void cppgpost2l(float part[], float q[], int npp, int noff, float qm,
                int idimp, int npmax, int nxv, int nypmx);

void cppgjpost2l(float part[], float cu[], float edges[], int npp,
                 int noff, int ihole[], float qm, float dt, int nx,
                 int ny, int idimp, int npmax, int nxv, int nypmx,
                 int idps, int ntmax, int ipbc);

void cppgrjpost2l(float part[], float cu[], float edges[], int npp,
                  int noff, int ihole[], float qm, float dt, float ci,
                  int nx, int ny, int idimp, int npmax, int nxv,
                  int nypmx, int idps, int ntmax, int ipbc);

void cppdsortp2yl(float parta[], float partb[], int npic[], int npp,
                  int noff, int nyp, int idimp, int npmax, int nypm1);

void cppcguard2xl(float fxy[], int myp, int nx, int ndim, int nxe,
                  int nypmx);

void cppaguard2xl(float q[], int myp, int nx, int nxe, int nypmx);

void cppacguard2xl(float cu[], int myp, int nx, int ndim, int nxe,
                   int nypmx);

void cppois23(std::complex<float> q[], std::complex<float> fxy[], 
			  int isign, std::complex<float> ffc[], float ax, 
			  float ay, float affp, float *we, int nx, int ny, 
			  int kstrt, int nyv, int kxp, int nyhd);

void cppcuperp2(std::complex<float> cu[], int nx, int ny, int kstrt, 
				int nyv, int kxp);

void cippbpoisp23(std::complex<float> cu[], std::complex<float> bxy[],
                  std::complex<float> ffc[], float ci, float *wm, int nx,
                  int ny, int kstrt, int nyv, int kxp, int nyhd);

void cppmaxwel2(std::complex<float> exy[], std::complex<float> bxy[],
                std::complex<float> cu[], std::complex<float> ffc[], float affp,
                float ci, float dt, float *wf, float *wm, int nx,
                int ny, int kstrt, int nyv, int kxp, int nyhd);

void cppemfield2(std::complex<float> fxy[], std::complex<float> exy[],
                 std::complex<float> ffc[], int isign, int nx, int ny,
                 int kstrt, int nyv, int kxp, int nyhd);

void cwpfft2rinit(int mixup[], std::complex<float> sct[], int indx, int indy,
                  int nxhyd, int nxyhd);

void cppfft2rxx(std::complex<float> f[], int isign, int mixup[],
                std::complex<float> sct[], int indx, int indy, int kstrt,
                int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                int nxyhd);

void cppfft2rxy(std::complex<float> g[], int isign, int mixup[],
                std::complex<float> sct[], int indx, int indy, int kstrt,
                int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                int nxyhd);


void cppfft2r3xx(std::complex<float> f[], int isign, int mixup[],
                 std::complex<float> sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd);

void cppfft2r3xy(std::complex<float> g[], int isign, int mixup[],
                 std::complex<float> sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd);

void cwppfft2r(std::complex<float> f[], std::complex<float> g[], 
			   std::complex<float> bs[], std::complex<float> br[],
			   int isign, int ntpose, int mixup[],
               std::complex<float> sct[], float *ttp, int indx, int indy,
               int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
               int kypd, int nxhyd, int nxyhd);

void cwppfft2r3(std::complex<float> f[], std::complex<float> g[], 
				std::complex<float> bs[], std::complex<float> br[], 
				int isign, int ntpose, int mixup[],
				std::complex<float> sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd);
