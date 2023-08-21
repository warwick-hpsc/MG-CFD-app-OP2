/* C header file for gpupfft2.cu */

void gpupfft2rrcuinit(int nx, int kypp, int ndim);

void gpupfft2cuinit(int kxpp, int ny, int ndim);

void gpupfft2rrcudel();

void gpupfft2cudel();

void gpupfft2rrcux(std::complex<float> f[], std::complex<float> bsm[], int isign,
                   int indx, int indy, int kstrt, int nvp, int kxp1,
                   int kyp, int nxh1d, int kypd);

void gpupfft2rrcuy(std::complex<float> g[], std::complex<float> brm[], int isign,
                   int indx, int indy, int kstrt, int nvp, int kxp1,
                   int kyp, int nyd);

void gpupfft2rrcuxn(std::complex<float> fn[], std::complex<float> bsm[], int isign,
                    int indx, int indy, int ndim, int kstrt, int nvp,
                    int kxp1, int kyp, int nxh1d, int kypd);

void gpupfft2rrcuyn(std::complex<float> gn[], std::complex<float> brm[], int isign,
                    int indx, int indy, int ndim, int kstrt, int nvp,
                    int kxp1, int kyp, int nyd);

void cgpuppsltpose(std::complex<float> f[], std::complex<float> g[], float ani,
                   int nx, int ny, int kxp, int kyp, int kstrt, int nxv,
                   int nyv);

void cgpuppsltposen(std::complex<float> fn[], std::complex<float> gn[], float ani,
                    int nx, int ny, int kxp, int kyp, int kstrt,
                    int ndim, int nxv, int nyv);
