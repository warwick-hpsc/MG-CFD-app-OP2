/* header file for gpplib2.c */

void cgppcacguard2l(std::complex<float> g_cu[], float g_cue[], float g_scs[],
                    float scs[], float scr[], int nx, int nyp,
                    int kstrt, int nvp, int ndim, int nxe, int nypmx,
                    int nxvh, int kypd);

void cgppcaguard2l(std::complex<float> g_q[], float g_qe[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int nxe, int nypmx, int nxvh, int kypd);

void cgppcbguard2l(std::complex<float> g_fxyz[], float g_fxyze[],
                   float g_scs[], float scs[], float scr[], int nx,
                   int nyp, int kstrt, int nvp, int ndim, int nxe,
                   int nypmx, int nxvh, int kypd);

void cwappfft2rcs(std::complex<float> g_f[], std::complex<float> g_g[],
                  std::complex<float> g_bsm[], std::complex<float> g_brm[],
                  std::complex<float> bsm[], std::complex<float> brm[], int isign,
                  int g_mixup[], std::complex<float> g_sct[], float ttp[],
                  int indx, int indy, int kstrt, int nvp, int kxpd,
                  int kyp, int nxhd, int nyd, int kypd, int nxhyd,
                  int nxyhd);

void cwappfft2rcsn(std::complex<float> g_fn[], std::complex<float> g_gn[],
                   std::complex<float> g_bsm[], std::complex<float> g_brm[],
                   std::complex<float> bsm[], std::complex<float> brm[], int isign,
                   int g_mixup[], std::complex<float> g_sct[], float ttp[],
                   int indx, int indy, int kstrt, int nvp, int ndim,
                   int kxpd, int kyp, int nxhd, int nyd, int kypd,
                   int nxhyd, int nxyhd);

void gpuppfft2rrcu(std::complex<float> g_f[], std::complex<float> g_g[],
                   std::complex<float> g_bsm[], std::complex<float> g_brm[],
                   std::complex<float> bsm[], std::complex<float> brm[], int isign,
                   float ttp[], int indx, int indy, int kstrt, int nvp,
                   int kxpd, int kyp, int nxhd, int nyd, int kypd);

void gpuppfft2rrcun(std::complex<float> g_fn[], std::complex<float> g_gn[],
                    std::complex<float> g_bsm[], std::complex<float> g_brm[],
                    std::complex<float> bsm[], std::complex<float> brm[], int isign,
                    float ttp[], int indx, int indy, int kstrt, int nvp,
                    int ndim, int kxpd, int kyp, int nxhd, int nyd,
                    int kypd);

void cgpporder2l(float g_ppart[], float g_ppbuff[], float g_sbufl[],
                 float g_sbufr[], int g_kpic[], int g_ncl[],
                 int g_ihole[], int g_ncll[], int g_nclr[], 
                 float sbufl[], float sbufr[], float rbufl[], 
                 float rbufr[], int ncll[], int nclr[], int mcll[],
                 int mclr[], float ttp[], int noff, int nyp, int kstrt,
                 int nvp, int idimp, int nppmx, int nx, int ny, int mx,
                 int my, int mx1, int myp1, int npbmx, int ntmax,
                 int nbmax, int *g_irc);

void cgpptpose(std::complex<float> g_bsm[], std::complex<float> g_btm[], 
               std::complex<float> sm[], std::complex<float> tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp);

void cgpptposen(std::complex<float> g_bsm[], std::complex<float> g_btm[],
                std::complex<float> sm[], std::complex<float> tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim);
