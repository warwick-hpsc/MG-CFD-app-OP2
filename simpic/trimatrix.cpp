#include "trimatrix.hpp"

TriMatrix::TriMatrix(int i)
{
  n = i;
  nc = n-1;
  tri_a = new Scalar[n];
  tri_b = new Scalar[n];
  tri_c = new Scalar[n];
  gam = new Scalar[n];

  //  new (tri_a) Scalar(0.0);

  for(i=0; i < n; i++)
    {
      tri_a[i] = tri_b[i] = tri_c[i] = 0.0;
    }

}

TriMatrix::~TriMatrix()
{
  delete [] tri_a;
  delete [] tri_b;
  delete [] tri_c;
  delete [] gam;
}

void TriMatrix::solve(Scalar *rhs, Scalar *ans, int nstrt)
{
  register int j;
  Scalar bet = tri_b[nstrt];

  ans[nstrt] = rhs[nstrt]/bet;
  for(j=nstrt+1; j<n; j++) {
    gam[j] = tri_c[j-1]/bet;
    bet = tri_b[j] - tri_a[j]*gam[j];
    ans[j] = (rhs[j] - tri_a[j]*ans[j-1])/bet;
  }

  for(j=nc-1; j>= nstrt; j--) ans[j] -= gam[j+1]*ans[j+1];
}

void TriMatrix::operator *= (Scalar f)
{
  for(int i=0; i < n; i++)
    {
      tri_a[i] *= f;
      tri_b[i] *= f;
      tri_c[i] *= f;
    }
}

void TriMatrix::print(FILE *fp)
{
  fprintf(fp, "\nMatrix: %d elements", n);
  for(int i=0; i < n; i++)
    {
      fprintf(fp, "\ni=%d; a=%g; b=%g; c=%g", 
	      i, tri_a[i], tri_b[i], tri_c[i]);
    }
  fprintf(fp, "\n");
}
