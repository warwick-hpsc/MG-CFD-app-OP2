#ifndef TRIMATRIX
#define TRIMATRIX

#include "functions.hpp"
#include "message.hpp"

// ITMAX : maximum number of iterations for nonlinearsolve (default)
#define ITMAX 1000

class TriMatrix: public Message  // tridiagonal matrix class
{
 private:
  int n, nc;  // n=# of rows; nc=n-1
  Scalar *gam;

  // matrix coefficients
  Scalar *tri_a, *tri_b, *tri_c;

  // variables for nonlinear solver
  Scalar *tri_b_copy, *rhs_copy;
  bool *fixed; // for Dirichlet boundary conditions
  int Lnorm;   // norm for evaluation of residual error
  bool nonlinearmode;

 public:
  TriMatrix(int i);
  ~TriMatrix();

  // functions for memory allocation/deallocation
  void alloc(int i);  // allocate memory for i-row matrix
  void dealloc(void); // deallocate matrix memory

  // accessor functions
  inline int get_n(void) { return n; }
  inline int get_nc(void) { return nc; }
  inline Scalar & a(int i) { return tri_a[i]; }
  inline Scalar & b(int i) { return tri_b[i]; }
  inline Scalar & c(int i) { return tri_c[i]; }
  inline bool is_fixed(int i) { return fixed[i]; }
  inline void set_fixed(int i, bool a=true) { fixed[i] = a; } 
  inline bool is_nonlinear(void) { return nonlinearmode; }
  
  // access matrix arrays -- avoid when possible
  inline Scalar * get_tri_a(void) { return tri_a; }
  inline Scalar * get_tri_c(void) { return tri_c; }
  inline Scalar * get_tri_b(void)
  {
    if(nonlinearmode) 
      {
	return tri_b_copy;
      }
    return tri_b;
  }

  // matrix operators
  TriMatrix& operator = (TriMatrix * RHS);
  void operator *= (Scalar f);
  void operator += (TriMatrix & RHS);

  // printing functions
  void print(FILE *fp=stdout);
  void print(char *fname);

  // matrix-vector multiply
  void multiply(Scalar *vec, Scalar *ans);

  // set matrix to zero
  void zero(void);

  // linear solve
  void solve(Scalar *rhs, Scalar *ans, int nstrt=0);

  // nonlinearsolve -- Newton-Raphson
  // (see Cartwright et al., Phys. Plasmas, Aug 2000)
  template <typename T>
  inline void nonlinearsolve(Scalar *rhs, Scalar *ans, 
			     Scalar (T::*f)(Scalar), 
			     Scalar (T::*fprime)(Scalar),
			     T * evaluator,
			     Scalar errmin = -1.,
			     int itmax=ITMAX,
			     Scalar (T::*error)(Scalar *rhs_new, 
						Scalar *ans)=0)



  {
    nonlinearmode = true;  // enter nonlinearmode

    // mindgame: to avoid compilation error
    Scalar err = 0;
    Scalar *btmp;
    
    // allocate copy of diagonal (if null)
    if(!tri_b_copy)
      {	
	tri_b_copy = new Scalar[n];
	rhs_copy = new Scalar[n];
      }

    // swap pointers
    btmp = tri_b;
    tri_b = tri_b_copy;
    tri_b_copy = btmp;

    // iterate to convergence
    int i;
    register int j;
    for(i=0; i < itmax; i++)
      {
	// copy variables
	memcpy(tri_b, tri_b_copy, n*sizeof(Scalar));
	memcpy(rhs_copy, rhs, n*sizeof(Scalar));

	// modify matrix coefficients
	for(j=0; j < n; j++)
	  {
	    if(!fixed[j])
	      {
		tri_b[j] -= (evaluator->*fprime)(ans[j]);
		rhs_copy[j] += (evaluator->*f)(ans[j])
		      - (evaluator->*fprime)(ans[j])*ans[j];	
	      }
	  }

	// trying to find out why n0const works and quasineutral doesn't,
	// JH, July 6, 2006	
#ifdef DEBUG
	printarray("tri_rhs.txt", n, rhs_copy);
#endif

	// solve linear system
	solve(rhs_copy, ans);
	
	// compute "true" rhs
	for(j=0; j < n; j++)
	  {
	    rhs_copy[j] = rhs[j];

	    rhs_copy[j] += (evaluator->*f)(ans[j]);
	    
	  }

	// calculate error
	
	if(evaluator && error) // evaluate error using S::*error
	  {
	    err = ((*evaluator).*error)(ans, rhs_copy);	   
	  }
	else // use built-in error evalution
	  {
	    Scalar err_loc;
	    Scalar rhs_loc;
	    Scalar denom;
	    for(j=1, err=0., denom=0.; j < nc; j++)
	      {
		err_loc = tri_a[j]*ans[j-1] 
		  + tri_b_copy[j]*ans[j] + tri_c[j]*ans[j+1];
		
		rhs_loc = rhs[j];
		
		rhs_loc += (evaluator->*f)(ans[j]);

		err_loc -= rhs_loc;
		
		err_loc = fabs(ppow(err_loc, Lnorm));
		denom += fabs(ppow(rhs_loc, Lnorm));
		
		err += err_loc;
	      }
	    
	    // normalize error
	    err = pow(double(err), double(1./Scalar(Lnorm))); 
	    denom = pow(double(denom), double(1./Scalar(Lnorm)));
	    err /= denom;	    
	  }

// trying to debug solver error tolerence, JH, Aug. 8, 2005
#ifdef DEBUG
	    fprintf(stderr, "i=%d,err=%g(%g)\n", 
		    i, err, errmin);
#endif

	    if(err < errmin) { break; }  // converged
      }    

    if(i == itmax)
      {
	// maximum iteration exceeded
	fprintf(stderr, "\nWARNING TriMatrix::nonlinearsolve() --\n");
	fprintf(stderr, "  %d iterations exceeded", itmax);
	fprintf(stderr, ", Residual error: %g ( > %g )\n", err, errmin);
      }
    // print # of iterations to convergence, JH, Aug. 8, 2005
#ifdef DEBUG
    else
      {
	fprintf(stderr, "\nTriMatrix::nonlinearsolve() -- ");
	fprintf(stderr, "converged in %d iterations (error=%g < %g)\n",
		i, err, errmin);
      }
#endif


    // swap pointers back
    btmp = tri_b;
    tri_b = tri_b_copy;
    tri_b_copy = btmp;

    // calculate true rhs??? NOT DONE?  JH, Aug. 11, 2005
    for(j=0; j < n; j++)
      {

      }
    nonlinearmode = false; // exit nonlinear mode
  }
};

#endif
