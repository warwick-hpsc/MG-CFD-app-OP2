void printarrays(const char *fname, int length, int narrays, Scalar **thearrays)
{
  int i, j;
  FILE *fp;

  fp = fopen(fname, "w");

  for(i=0; i < length; i++)
    {
      for(j=0; j < narrays; j++)
	{
	  fprintf(fp, "%g ", thearrays[j][i]);
	}
      fprintf(fp, "\n");
    }

  fclose(fp);
}

void printarrays(const char *fname, int length, int narrays, ...)
{
  va_list arguments;        //A place to store the list of arguments
  Scalar **thearrays;

  va_start(arguments, narrays);

  thearrays = new Scalar*[narrays];

  for(int i=0; i<narrays; i++)   //Loop until all numbers are added
    {
      thearrays[i] = va_arg(arguments, Scalar *);
    }

  printarrays(fname, length, narrays, thearrays);

  va_end(arguments);      //Cleans up the list
  delete [] thearrays;  
}

void printtarray(int ndiag, int arraysize, Scalar incr, int close=0)
{
  register int i;

  static int init = 0;
  static FILE *fp[NDIAGNOSTICS];
  
  char ffname[STR_SIZE];

  if(!init)
    {

      #ifdef DEBUG
      //      fprintf(stderr, "\ninitializing printtarray():\n");
      #endif 

      for(i=0; i < NDIAGNOSTICS; i++)
	{

	  #ifdef DEBUG
	  //	  fprintf(stderr, "\n!!!!diagnostic %d:\n", i);
	  #endif

	  sprintf(ffname, "%s%s", diags[i], fext);

	  #ifdef DEBUG
	  //	  fprintf(stderr, "\n%s, %s\n", diags[i], fext);
	  //	  fprintf(stderr, "\n%s\n", ffname);
	  #endif

	  fp[i] = fopen(ffname, "w");

	  #ifdef DEBUG
	  //	  fprintf(stderr, "diagnostic %d of %d\n", i, NDIAGNOSTICS);
	  //	  fprintf(stderr, "-- %s\n", ffname);
	  //	  fprintf(stderr, "\n...done with this diagnostic!\n");
	  #endif
	}

      init = 1;
    }

  if(close)
    {
      for(i=0; i < NDIAGNOSTICS; i++)
	{
	  fclose(fp[i]);
	}
      return;
    }


  for(i=0; i < arraysize; i++)
    {
      fprintf(fp[ndiag], "%g %g %g\n", tt, incr*i, thediags[ndiag][i]);
    }
  fprintf(fp[ndiag], "\n");

}

void printphase(char *fname)
{
  register int i;
  register int j;

  FILE *fp;

  fp = fopen(fname, "w+");

  for(i=0, j=0; i < npart; i++, j+=2)
    {
      fprintf(fp, "%g %g %g\n", tt, pdata[j], pdata[j+1]);
    }
  fprintf(fp, "\n");

  fclose(fp);
}

void printt(int close=0)
{
  static int init = 0;
  static FILE *fp;
  char fname[STR_SIZE];
  
  if(!init)
    {
      #ifdef DEBUG
      //      fprintf(stderr, "\nInitializing printt...\n");
      #endif

      sprintf(fname, "nt%s.dat", procname);

      fp = fopen(fname, "w");


      init = 1;
    }

 if(close)
    {
      fclose(fp);
      return;
    }

  fprintf(fp, "%g %d\n", tt, npart);
}
 
