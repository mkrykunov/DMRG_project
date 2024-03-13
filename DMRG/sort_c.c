void sort_c(int n, int SMALL, double d[], double **v)
{
	int i, j, k;
	double p;

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
	      for (j = 0; j < n; j++)
	      {
		 p = v[j][i];
		 v[j][i] = v[j][k];
		 v[j][k] = p;
	      }
	   }
	}
}
