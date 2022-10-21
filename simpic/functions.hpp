#ifndef	FUNCTIONS_H
#define	FUNCTIONS_H

#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <cmath>
#include <vector>

// mindgame: this file can be generalized into a kind of math library.

//#define Scalar float 
#define EPS 3.0e-7
#define ITMAX 50
#define STR_SIZE 501

// definitions
#define sqr(x) ((x)*(x))
// added by JK, 2017-09-22
#ifndef XG_SCALAR_DOUBLE
#define pow(x,y) pow((double) (x), (double) (y))
#endif

// forward declarations

const long frandseed=31207321;

inline Scalar gammln(Scalar xx)
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}


//	return a random number between 0 and 1
inline float frand2() {float f = rand(); return f/((float)RAND_MAX+1.0);}

/************************************************************************/
/* frand() returns values 0 through 1.                                  */
/* From "Random number generators: good ones are hard to find", S. Park */
/* and K. Miller, Communications of ACM, October 1988, pp 1192-1201.    */
/* This is from page 1195, and is to work on any system for which       */
/* maxint is 2**31-1 or larger. Due earlier to Schrage, as cited by P&M.*/
/*                                                                      */
/* Note: OK values for iseed are 1 through 2**31-2. Give it 0 or 2*31-1 */
/* and it will return the same values thereafter!                       */
/*                                                                      */
/* C version 6/91, Bruce Langdon.                                       */
/*                                                                      */
/* Algorithm replaces seed by mod(a*seed,m). First represent            */
/* seed = q*hi + lo.  Then                                              */
/* a*seed = a*q*hi + a*lo = (m - r)*hi + a*lo = (a*lo - r*hi) + m*hi,     */
/* and new seed = a*lo - r*hi unless negative; if so, then add m.       */

// extern long frandseed;

inline float frand()
{
  long a = 16807, m = 2147483647, q = 127773, r = 2836;
  long hi, lo;
  float fnumb;
  static long frandseed=31207321; 

	// moved seed to frandseed
  //  fnumb = 2;
  // the random number returned is 0=<fnumb<1
  do {
	 hi = frandseed/q;
	 lo = frandseed - q*hi;
	 frandseed = a*lo - r*hi;
	 // "seed" will always be a legal integer of 32 bits (including sign)
	 // if(seed <= 0) seed = seed + m;
	 if(frandseed < 0) frandseed = frandseed + m;
	 fnumb = frandseed/(2147483647.0);
  } while (fnumb>=1);

  return(fnumb);
}

inline float normal()
{
  // returns a normally distributed deviate with zero mean and unit variance
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	
	if(iset==0){
		do{
			v1=2.0*frand()-1.0;
			v2=2.0*frand()-1.0;
			r = v1*v1+v2*v2;
		} while (r>=1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset =1;
		r=v2*fac;
	}
	else{
		iset=0;
		r = gset;
	}
	return(r);
}

/*returns only a positive one in normal distribution.*/
inline float half_normal()
{
// mindgame
//	return (sqrt(2.0*fabs(log(fabs(frand()+1e-7)))));
  Scalar R = 1-frand();
  Scalar foo = (sqrt(-2.0*log(R)));
  return foo;
}


//	return a randomly signed 1
inline int rsign() {int i = rand()%2; if (i == 0) return -1; else return 1;}



// return a bit reversed number
inline float revers(int num, int n)
{
  float power, rev;
  int inum, iquot, irem;
  
  rev = 0.;
  inum = num;
  power = 1.;
  
  do
	 {
		iquot = inum/n;
		irem = inum - n*iquot;
		power /= n;
		rev += irem*power;
		inum = iquot;
	 } 
  while (inum > 0);
  
  return (rev);
}



/* Returns bit reversed num in base 2 */

inline float revers2(unsigned int num)
{
  double f=0.5, sum=0.0;
  
  while(num)
		{
			if (num & 1) sum += f;     /* is 1st bit set? */
			f *= 0.5;
			num >>= 1;		       /* fast divide by 2 */
		}
  return (sum);
}

inline float base2()
{ 
  static int counter2=0;
  return(revers(counter2++, 2));
}

inline float revers_3()
{ 
  static int iset=0;
  static int counter=0;
  if(iset==0){
    iset=1;
    return(revers(counter++, 3));
  }
  else if(iset==1){
    iset=2;
    return(revers(counter, 5));
  }
  else{
    iset=0;
    return(revers(counter, 7));
  } 
}

inline float revers_normal()
{
  /*returns a normally distributed deviate with zero mean and unit variance.*/
  static int iset=0;
  static int count=0;
  static float gset;
  float fac,r,v1,v2;
	
  if(iset==0){
    do{
      v1=2.0*revers(count++,2)-1.0;
      v2=2.0*revers(count,5)-1.0;
      r = v1*v1+v2*v2;
    } while (r>=1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset =1;
    r=v2*fac;
  }
  else{
    iset=0;
    r = gset;
  }
  return(r);
}

// mindgame
/*returns only a positive one in normal distribution.*/
inline float half_revers_normal()
{
// mindgame
//	return (sqrt(2.0*fabs(log(fabs(frand()+1e-7)))));
	return (sqrt(-2.0*log(1-base2())));
}


inline Scalar rtbis(double func(double), Scalar value, int moments,
		    Scalar x1, Scalar x2)
{
  int j;
  Scalar dx,f,fmid,xmid,rtb;

  f=func(x1)-value;
  fmid=func(x2)-value;
  if (f*fmid > 0.0) fprintf(stderr, "\nrtbis: Root must be bracketed.");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=ITMAX;j++) {
    fmid=func(xmid=rtb+(dx *= 0.5))-value;
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < EPS || fmid == 0.0) {
      return rtb;
    }
    f=fmid;
  }
  fprintf(stderr, "\nrtbis: Maximum number of iterations exceeded %e.", 
	  value/func(x2));
  return 0.0;                     /* Never get here. */
}

#undef ITMAX
#undef EPS

//------------------------------------------------------------------
// Support for error functions and complementary error functions.
//	From Numerical Recipes in C... These have been ANSI-fied for
//	use in c++.

#define ITMAX 100
#define EPS 3.0e-7
inline void gcf(Scalar* gammcf, Scalar a, Scalar x, Scalar* gln)
{
	int n;
	Scalar gold=0.0,g,fac=1.0,b1=1.0;
	Scalar b0=0.0,anf,ana,an,a1,a0=1.0;
	//	Scalar gammln(Scalar x);
//	void nrerror();

	*gln=gammln(a);
	a1=x;
	for (n=1;n<=ITMAX;n++) {
		an=static_cast<Scalar>(n);
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < EPS) {
				*gammcf=exp(-x+a*log(x)-(*gln))*g;
				return;
			}
			gold=g;
		}
	}
//  mindgame: to avoid compilation warning
	*gammcf=0;
//	nrerror("a too large, ITMAX too small in routine GCF");
}
#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7
inline void gser(Scalar* gamser,Scalar a, Scalar x, Scalar* gln)
{
	int n;
	Scalar sum,del,ap;
	Scalar gammln(Scalar x);
//	void nrerror();

	*gln=gammln(a);
	if (x <= 0.0) {
//		if (x < 0.0) nrerror("x less than 0 in routine GSER");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
//		nrerror("a too large, ITMAX too small in routine GSER");
//  mindgame: to avoid compilation warning
		*gamser=0.0;
		return;
	}
}
#undef ITMAX
#undef EPS

inline Scalar gammp(Scalar a, Scalar x)
{
//	Scalar gamser,gammcf,gln;
	Scalar gam,gln;
	//	void gser(Scalar* gamser,Scalar a, Scalar x, Scalar* gln);
	//	void gcf(Scalar* gammcf, Scalar a, Scalar x, Scalar* gln);
//	void nrerror();

//	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMP");
	if (x < (a+1.0)) {
//		gser(&gamser,a,x,&gln);
//		return gamser;
		gser(&gam,a,x,&gln);
		return gam;
	} else {
//		gcf(&gammcf,a,x,&gln);
//		return (1.0-gammcf);
		gcf(&gam,a,x,&gln);
		return (1.0-gam);
	}
}

inline Scalar gammq(Scalar a, Scalar x)
{
	Scalar gamser,gammcf,gln;
//	void gcf(),gser(),nrerror();

//	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine GAMMQ");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return (1.0-gamser);
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

inline Scalar erf(Scalar x)
{
//	Scalar gammp();

	return (x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x));
}

inline Scalar erfc(Scalar x)
{
//	Scalar gammp(),gammq();

	return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}

inline Scalar erfcc(Scalar x)
{
	Scalar t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return  x >= 0.0 ? ans : 2.0-ans;
}

inline Scalar bessj0(Scalar x)
{
	Scalar ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
				 +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
				 +y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
           +y*(-0.6911147651e-5+y*(0.7621095161e-6
           -y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

inline Scalar bessj1(Scalar x)
{
	Scalar ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
				 +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

inline float normal2()
{
  // returns a normally distributed deviate with 
  // zero mean and 1/Sqrt(2) variance.
  // The way that we have defined v thermal 
  // this returns gaussian with a v thermal
	/* of one. */
	static int iset=0;
	static float gset;
	float fac,r2,v1,v2;
	
	if(iset==0){
		do{
			v1=2.0*frand()-1.0;
			v2=2.0*frand()-1.0;
			r2 = v1*v1+v2*v2;
		} while (r2>=1.0);
		fac=sqrt(-log(r2)/r2);
		gset=v1*fac;
		iset =1;
		r2=v2*fac;
	}
	else{
		iset=0;
		r2 = gset;
	}
	return(r2);
}
/*
float normal()
{
	returns a normally distributed deviate with zero mean and unit variance.
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	
	if(iset==0){
		do{
			v1=2.0*frand()-1.0;
			v2=2.0*frand()-1.0;
			r = v1*v1+v2*v2;
		} while (r>=1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset =1;
		r=v2*fac;
	}
	else{
		iset=0;
		r = gset;
	}
	return(r);
}
*/

// returns the positive integer power of a Scalar
inline Scalar ppow(Scalar a, int pow)
{
  Scalar val = 1.;
  for(int i=0; i < pow; i++)
    {
      val *= a;
    }
  return val;
}

inline Scalar ipow(Scalar a, int pow)
{
  bool negflag = 0;

  if(pow < 0) { negflag = 1; pow = -pow; }

  Scalar val = ppow(a, pow);

  if(negflag) { val = 1./val; };
  return val;
}

#ifndef _STANDALONE_LOADLIB

/*  This function takes two ordered pairs of integers,
which represent a line, and
returns an array of ordered pairs of integers which 
represent horizontal or vertical line segments which
will interpolate the given line onto a mesh. The array
ends when four 0 integers are given. */

inline int *segmentate(int j1, int k1, int j2, int k2) 
{
  int j;
  int jl,kl,jh,kh;
  int *segments;
  Scalar m,y;
  int kt;
  int j1s, k1s, j2s,k2s;  /*  the line segment coordinates */
  int segcount = 0;  /* number of segments. */
  /* reorder the input such that jl,kl is the leftmost point */
  if(j1<j2) { 
	  jl=j1;kl=k1;
	  jh=j2;kh=k2;
  }
  else {
	  jl=j2;kl=k2;
	  jh=j1;kh=k1;
  }

  m = (kh -kl)/static_cast<Scalar>(jh - jl);
  kt = kl;

  segments=new int[4*(3+ (jh-jl))];
  memset(segments,0,4*(3+(jh-jl))*sizeof(int));

  /* make sure the first point is in the array */
  segments[0]=jl;
  segments[1]=kl;

  /* get all the horizontal segments. */
  for(j=jl;j<jh;j++) {
	  y = m * (j - jl + 0.5) + kt;
	  k1s = (int) (y + 0.5);
	  k2s = k1s;
	  j1s = j;
	  j2s = j+1;
	  segments[4*segcount+2]=j1s;
	  segments[4*segcount+3]=k1s;
	  segments[4*segcount+4]=j2s;
	  segments[4*segcount+5]=k2s;
	  segcount++;
  }

  /* make sure the last point is in the array */
  segments[4*segcount+2]=jh;
  segments[4*segcount+3]=kh;
  return segments;  // make sure to delete this later
}

inline Scalar LineDist(Scalar A1, Scalar A2, 
		       Scalar B1, Scalar B2, Scalar tA1, Scalar tA2) {
	Scalar dy = A2 - B2;
	Scalar dx = A1 - B1;
	return fabs( dx * tA2 - dy * tA1 +dy * A1 - dx * A2 )/ sqrt( dx*dx + dy*dy);
}  

inline int nlines(char *filename)
{
  int i=0;
  char buffer[STR_SIZE];
  FILE *fp;

  if((fp = fopen(filename, "r")) == NULL)
    {
      return -1;
    }

  while(fgets(buffer, STR_SIZE, fp))
    {
      i++;
    }
  fclose(fp);
  return i;
}

// returns an array of Scalars that are equally spaced on a log plot
// -- must be deleted by caller (using delete)
inline Scalar * logspacing(Scalar Emin, Scalar Emax, int n)
{
  Scalar *a;
  a = new Scalar[n+1];
  Scalar C = log(pow(Emax/Emin, Scalar(1./Scalar(n))));
  a[0] = Emin;
  for(int i=0; i < n; i++)
    {
      a[i+1] = a[i]*exp(C);
    }
  return a;
}

// returns an array of Scalars that are equally spaced on a linear plot
// -- must be deleted by caller (using delete)
inline Scalar * linearspacing(Scalar Emin, Scalar Emax, int n)
{
  Scalar *a;
  a = new Scalar[n+1];
  Scalar dx = (Emax-Emin)/static_cast<Scalar>(n);
  a[0] = Emin;
  for(int i=0; i < n; i++)
    {
      a[i+1] = a[i] + dx;
    }

  return a;
}

// mindgame: linear interpolation fn.
//           assuming p1.first != p0.first
template <typename Type1, typename Type2>
inline Type2 linear_interpolation(std::pair<Type1, Type2> const &p0,
				  std::pair<Type1, Type2> const &p1,
				  Type1 const &x)
{
  return (p1.second-p0.second)*static_cast<Type2>(x-p0.first)
	 /static_cast<Type2>(p1.first-p0.first) + p0.second; 
}

// mindgame: linear interpolation fn.
template <typename Type1, typename Type2>
inline Type2 safe_linear_interpolation(std::pair<Type1, Type2> const &p0,
				       std::pair<Type1, Type2> const &p1,
				       Type1 const &x)
{
  if (static_cast<Type2>(p1.first-p0.first) != 0)
    linear_interpolation(p0, p1, x);
  else return static_cast<Type2>(0);
}

void obsolete_printarrays(const char *name, int arrsize, int narrays, ...)
{
  int i;
  unsigned int j;
  std::vector<Scalar *> thearrays;
  va_list argp;

  FILE * thefp = fopen(name, "w");
  if(!thefp) 
    {
      fprintf(stderr, "File '%s' could not be opened", name);
      return;
    }

  va_start(argp, narrays);
  for(i=0; i < narrays; i++)
    {
      thearrays.push_back(va_arg(argp, Scalar *));
    }

  for(i=0; i < arrsize; i++)
    {
      for(j=0; j < thearrays.size(); j++)
	{
	  fprintf(thefp, "%g ", thearrays[j][i]);
	}
      fprintf(thefp, "\n");
    }

  fclose(thefp);
}

// printarray: front-end to printarrays for a single array
inline void printarray(const char *name, int arrsize, Scalar *arr)
{
  printarrays(name, arrsize, 1, arr);
}
#endif

#endif

