// fields
#include <algorithm> 
#define LAPLACE

TriMatrix *A;

Scalar rhsV(Scalar t)
{
  return 0.;
}

Scalar lhsV(Scalar t)
{
  return lhsvoltage;
}

void allocate_arrays(int dealloc=0)
{
  if(dealloc)
    {
      delete [] narray;
      delete [] phiarray;
      delete [] Earray;
      delete [] nback;

      return;
    }
  
  narray = new Scalar[ng];
  phiarray = new Scalar[ng];
  Earray = new Scalar[2*npart + 1];
  nback = new Scalar[ng];

  std::fill_n(Earray, 2*npart + 1, 0); 

  
}

void allocate_fieldbuffers(int dealloc=0)
{
  if(dealloc)
    {
      delete [] fsbuffer;
      delete [] frbuffer;
      return;
    }

  fsbuffer = new Scalar[2*nproc];
  frbuffer = new Scalar[2*nproc];
}

void sumLaplace(Scalar *pphh)
{
  register int i;
  Scalar rv, lv, frac, coord, xlocal;

  rv = rhsV(t);
  lv = lhsV(t);

  for(i=0, xlocal = xl; i < ng; i++, xlocal +=dx)
    {
      frac = xlocal/L;
      pphh[i] += frac*rv + (1. - frac)*lv;
    }
}

void setcoeffs(Scalar scale)
{
  A = new TriMatrix(ng);

  // lhs BC
  if(rank == 0)
    {
      A->a(0) = 0.0;
      A->b(0) = 1.0;
      A->c(0) = 0.0;
    }
  else
    {
      //      A->a(0) = ;
      A->b(0) = -scale*(1.0 + dx/xl);
      A->c(0) = scale;
       
    }
  
  // rhs BC
  if(rank == last)
    {
      A->a(nc) = 0.0;
      A->b(nc) = 1.0;
      A->c(nc) = 0.0;
    }
  else
    {
      A->a(nc) = scale;
      A->b(nc) = -scale*(1.0 + dx/(L - xr));
    }

  for(int i=1; i < nc; i++)
    {
      A->a(i) = scale;
      A->b(i) = -2.*scale;
      A->c(i) = scale;
    }

}

void init_fields()
{
  allocate_arrays();
  allocate_fieldbuffers();
  
  memset(nback, 0, ng*sizeof(Scalar));

  setcoeffs(-epsilon/(q*dx*dx));
}

void communicate_fields(void)
{
  int i, j, k;
  Scalar lcoord, rcoord;
  Scalar xval, phival;

  Scalar ttemp[2];

#ifdef DEBUG
  //  fprintf(stderr, "proc %d entering communicate_fields\n", rank);
#endif

  if(nproc > 1)
    {
      ttemp[0] = phiarray[0];
      ttemp[1] = phiarray[nc];

      MPI_Allgather(ttemp, 2, MPI_SCALAR, frbuffer, 2, MPI_SCALAR, custom_comm);


      for(i=0, j=0; i < nproc; i++, j+=2)
	{
#ifdef DEBUG
	  //	  fprintf(stderr, "rank %d %d: %g %g;\n", 
	  //		  rank, i, frbuffer[j], frbuffer[j+1]);
#endif
	}

      // compute contribution from lhs
      for(i=0, j=0; i < rank; i++, j+=2)
	{
	  phival = frbuffer[j+1];
	  lcoord = L*(i+1)/Scalar(nproc);
	  lcoord = L - lcoord;
	  for(k=0, xval=L-xl; k < ng; k++, xval -= dx)
	    {
	      phiarray[k] += phival*(xval/lcoord);
	    }
	}
      
      // compute contribution from rhs
      for(i=rank+1; i < nproc; i++)
	{
	  j = 2*i;
	  phival = frbuffer[j];
	  rcoord = L*i/Scalar(nproc);
	  for(k=0, xval=xl; k < ng; k++, xval += dx)
	    {
	      phiarray[k] += phival*(xval/rcoord);
	    }
	}

    }

#ifdef DEBUG
  //      fprintf(stderr, "proc %d done with communicate_fields\n", rank);
#endif
}

void advancefields(Scalar ddt)
{
  Scalar nlold, nrold;

  starttime(FIELDS);

  // modify density array
  if(rank == 0)
    {
      nlold = narray[0];
      narray[0] = 0.;
    }
  else
    {
      narray[0] *= 2;
    }
  if(rank == last)
    {
      nrold = narray[nc];
      narray[nc] = 0.;
    }
  else
    {
      narray[nc] *= 2;
    }

  A->solve(narray, phiarray);

  // restore density array
  if(rank == 0)
    {
      narray[0] = 2*nlold;
    }
  if(rank == last)
    {
      narray[nc] = 2*nrold;
    }
  
  endtime(FIELDS);
  starttime(MPITIME);
  communicate_fields();
  endtime(MPITIME);
  starttime(FIELDS);

  sumLaplace(phiarray);
  gradient(Earray, phiarray, ng, -0.5/dx);

  // "fix up" endpoints
  Earray[0] = -(phiarray[1] - phiarray[0])/dx;
  Earray[nc] = -(phiarray[nc] - phiarray[nc-1])/dx;
    

  endtime(FIELDS);
}
