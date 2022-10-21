// functions for timing the program

double ttemp[NTIMES];

inline void starttime(int type)
{
  ttemp[type] = MPI_Wtime();
}

inline void endtime(int type)
{
  tttimes[type] += (MPI_Wtime() - ttemp[type]);
}

// assumes FILE *fp is open already
void printtimes(FILE *fp)
{
  int i;
  Scalar sume, ttot;

  ttot = tend - tstart;

  fflush(fp);
  fprintf(fp, "\nTimings:\n");
  for(i=0, sume = 0.; i < NTIMES; i++)
    {
      fprintf(fp, "%s time: %g [%.1f%%]\n", names[i], tttimes[i],
	      100.*tttimes[i]/ttot);
      sume += tttimes[i];
    }
  fprintf(fp, "\nTotal time: %g (%g [%.1f%%])\n", ttot, sume, 100.*sume/ttot);
  fprintf(fp, "%g (%g) particles/second\n", nparttot/ttot, 
	  nparttot/tttimes[PARTICLES]);
  fflush(fp);
}
