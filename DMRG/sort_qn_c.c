void sort_qn_c(int n, int SMALL, double d[], double *v, double qn[], long index[], int index2[])
{
	int i, j, k;
	long ltmp;
	double p, tmp;

	for (i = 0; i < n-1; i++)
	{
	   k = i;
	   p = d[i];
	   for (j = i+1; j < n; j++)
	   {
	      if (SMALL)
	      {
		 if (d[j] <= p)
		 {
		    k = j;
		    p = d[j];
		 }
	      }
	      else
	      {
		 if (d[j] >= p)
		 {
		    k = j;
		    p = d[j];
		 }
	      }
	   }

	   if (k != i)
	   {
	      d[k] = d[i];
	      d[i] = p;

	      tmp = qn[i];
	      qn[i] = qn[k];
	      qn[k] = tmp;

              ltmp = index[i];
              index[i] = index[k];
              index[k] = ltmp;

	      index2[i] = k;

	      for (j = 0; j < n; j++)
	      {
                 p = v[i * n + j];
                 v[i * n + j] = v[k * n + j];
                 v[k * n + j] = p;
	      }
	   }
	}

        for (i = 0; i < n-1; i++)
        {
	   k = index2[i];

           if (k != i)
           {
              for (j = 0; j < n; j++)
              {
                 p = v[j * n + i];
                 v[j * n + i] = v[j * n + k];
                 v[j * n + k] = p;
              }
	   }
	}

}
