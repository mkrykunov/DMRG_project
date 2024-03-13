#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "Typedef.h"


#define FALSE	0
#define TRUE	1
#define sign(a, b)	( (b < 0) ? -fabs(a) : fabs(a) )


void Hx(SparseMatrix<T>& SA, T* HU, T* a, T* b)
{
	for ( int i = 0; i < SA.getSize(); i++ )
		a[i] = b[i] * SA[i];

	for ( int i = 0; i < SA.getNEL(); i++ )
	{
		a[SA.r(i)] += HU[i] * b[SA.c(i)];
		a[SA.c(i)] += HU[i] * b[SA.r(i)];
	}

	return;
}

void Davidson(int iEig, SparseMatrix<T>& SA, Vector<T>& En, Vector<T>& v, int iDim, T Eps)
{
	T	Zero = 1.0E-38;
	T	tmp, tmp1;
	int	nConv, iCD, iTop, iCount, Iter, itmp;
	int	i, j, k, lEnd, N = SA.getSize(), NE = SA.getNEL();
	T*	HU = SA.getVec() + N;
	char	str[40];

	iCount = Iter = nConv = 0;

	Matrix<T>	b(iDim,N), h(iDim,N), r(iEig,N), vtmp(N,iDim);
	Matrix<T>	a(iDim,iDim);
	Vector<T>	Ttmp(iDim*iDim); // , a(iDim*iDim);
	Vector<T>	Mu(iDim), ar(iDim), e(iDim);
	Vector<int>	IMinHD(iDim), lDisV(iDim);

	for ( i = 0; i < iEig; i++ )
		IMinHD[i] = 0;

	for ( i = iEig; i < iDim; i++ )
		ar[i] = 1.E3;

	for ( j = 0; j < N; j++ )
		h[0][j] = SA[j];

	for ( i = 0; i < iEig; i++ )
	{
		tmp = h[0][i];
		k = IMinHD[i] = i;
		for ( j = i+1; j < N; j++ )
		{
			if ( h[0][j] < tmp )
			{
				tmp = h[0][j];
				k = IMinHD[i] = j;
			}
		}
		if ( k != i ) 
		{
			h[0][k] = h[0][i];
			h[0][i] = tmp;
		}
	}

	for ( i = 0; i < iDim; i++ )
	{
		for ( j = 0; j < N; j++ )
		{
			b[i][j] = 0.0;
			h[i][j] = 0.0;
		}
	}

	for ( i = 0; i < iEig; i++ )
		b[i][IMinHD[i]] = 1.0;

	for ( i = 0; i < NE; i++ )
	{
		for ( j = 0; j < iEig; j++ )
		{
			if (SA.r(i) == IMinHD[j])
			{
				h[j][SA.c(i)] = HU[i];
			}
			else if (SA.c(i) == IMinHD[j])
			{
				h[j][SA.r(i)] = HU[i];
			}
		}
	}

	for ( i = 0; i < iEig; i++ )
		h[i][IMinHD[i]] = SA[IMinHD[i]];

	for ( i = 0; i < iEig; i++ )
	{
		for ( j = i; j < iEig; j++ )
		{
			Ttmp[i*iDim+j] = 0.0;
			for ( k = 0; k < N; k++ )
				Ttmp[i*iDim+j] += h[i][k] * b[j][k];
			Ttmp[j*iDim+i] = Ttmp[i*iDim+j];
		}
	}

	iCD = iEig;
	
	while (TRUE)
	{
		for ( i = 0; i < iCD; i++ )
		{
			for ( j = i; j < iCD; j++ )
			{
				a[i][j] = Ttmp[i*iDim+j];
				a[j][i] = Ttmp[j*iDim+i];
			}
		}

                tred2_c(a.getMat(), iCD, Mu.getVec(), e.getVec());
                tqli_c(Mu.getVec(), e.getVec(), iCD, a.getMat());

                int SMALL = 1;

                sort_c(iCD, SMALL, Mu.getVec(), a.getMat());

		Iter++;
		lEnd = TRUE;

		for ( i = 0; i < iEig; i++ )
		{
			for ( j = 0; j < N; j++ )
			{
				r[i][j] = 0.0;
				for ( k = 0; k < iCD; k++ )
					r[i][j] += (h[k][j] - Mu[i] * b[k][j]) * a[k][i];
			}
			ar[i] = 0.0;
			for ( j = 0; j < N; j++ )
				ar[i] += r[i][j] * r[i][j];
			ar[i] = sqrt(ar[i]);
			lEnd = ( lEnd && (ar[i] <= Eps) ) ? TRUE : FALSE;
		}

//		cout << "Number of products: " << iCount << endl;
//		_getche();

		double min = 0.0;

		if (iCount % 10 == 0)
		{
			min = ar[0];
			for (i = 0; i < iEig; i++)
				if (min > ar[i])	min = ar[i];

			sprintf(str, "Number of products: %d, min = %.3le", iCount, min);
			Trace(str);
		}

		if (lEnd)	break;

		iTop = (iCD+iEig > iDim) ? iEig : iCD;
		lDisV[0] = FALSE;

		for ( i = 0; i < iTop; i++ )
		{
			for ( j = 0; j < N; j++ )
			{
				vtmp[j][i] = 0.0;
				for ( k = 0; k < iCD; k++ )
					vtmp[j][i] += b[k][j] * a[k][i];
			}
		}

		for ( i = 0; i < iTop; i++ )
		{
			for ( j = 0; j < N; j++ )
				b[i][j] = vtmp[j][i];
		}

		for ( i = 0; i < iTop; i++ )
		{
			for ( j = 0; j < N; j++ )
			{
				vtmp[j][i] = 0.0;
				for ( k = 0; k < iCD; k++ )
					vtmp[j][i] += h[k][j] * a[k][i];
			}
		}

		for ( i = 0; i < iTop; i++ )
		{
			for ( j = 0; j < N; j++ )
				h[i][j] = vtmp[j][i];
		}

		lDisV[0] = TRUE;

		if ( (Iter % 15) == 0 )
		{
			for ( i = 0; i < iTop; i++ )
			{
				for ( j = i+1; j < iTop; j++ )
				{
					tmp = 0.0;
					for ( k = 0; k < N; k++ )
						tmp += b[i][k] * b[j][k];
					if ( fabs(tmp) > N*sqrt(Zero) )
					{
						lDisV[0] = TRUE;
						tmp1 = 0.0;
						if ( ar[i] <= ar[j] )
						{
							for ( k = 0; k < N; k++ )
							{
								b[j][k] -= tmp * b[i][k];
								tmp1 += b[j][k] * b[j][k];
								h[j][k] -= tmp * h[i][k];
							}
							if ( fabs( sqrt(tmp1)-1.0 ) > N*sqrt(Zero)*1.E-5 )
							{
								for ( k = 0; k < N; k++ )
								{
									b[j][k] /= sqrt(tmp1);
									h[j][k] /= sqrt(tmp1);
								}
							}
						}
						else
						{
							for ( k = 0; k < N; k++ )
							{
								b[i][k] -= tmp * b[j][k];
								tmp1 += b[i][k] * b[i][k];
								h[i][k] -= tmp * h[j][k];
							}
							if ( fabs( sqrt(tmp1)-1.0 ) > N*sqrt(Zero)*1.E-5 )
							{
								for ( k = 0; k < N; k++ )
								{
									b[i][k] /= sqrt(tmp1);
									h[i][k] /= sqrt(tmp1);
								}
							}
						}
					}
				}
			}
		}

		for ( i = 0; i < iTop; i++ )
		{
			for ( j = i; j < iTop; j++ )
			{
				Ttmp[i*iDim+j] = 0.0;
				for ( k = 0; k < N; k++ )
					Ttmp[i*iDim+j] += h[i][k] * b[j][k];
				Ttmp[j*iDim+i] = Ttmp[i*iDim+j];
			}
		}

		iCD = iTop;
		itmp = iCD;

		for ( i = 0; i < iEig; i++ )
		{
			if ( ar[i] <= Eps )	continue;
			for ( j = 0; j < N; j++ )
			{
				if ( fabs( SA[j] - Mu[i] ) > sqrt(Zero) )
					b[itmp][j] = -r[i][j] / ( SA[j] - Mu[i] );
				else
					b[itmp][j] = -r[i][j] / sign( sqrt( (T)Zero ), SA[j] - Mu[i] );
			}
			itmp++;
		}

		for ( i = iCD; i < itmp; i++ )
		{
			for ( k = 0; k < i; k++ )
			{
				tmp = 0.0;
				for ( j = 0; j < N; j++ )
					tmp += b[k][j] * b[i][j];
				tmp1 = 0.0;
				for ( j = 0; j < N; j++ )
				{
					b[i][j] -= tmp * b[k][j];
					if ( fabs(b[i][j]) <= sqrt(Zero) )	b[i][j] = 0.0;
					tmp1 += b[i][j] * b[i][j];
				}
				
				lDisV[i] = ( sqrt(tmp1) <= sqrt(Zero) ) ? TRUE : FALSE;
				if ( lDisV[i] )
				{
					for ( j = 0; j < N; j++ )
						b[i][j] = 0.0;
				}
				else
				{
					for ( j = 0; j < N; j++ )
					{
						b[i][j] = (fabs(b[i][j]) > Zero) ? b[i][j] / sqrt(tmp1) : 0.0;
					}
				}
			}
		}

		k = itmp;

		for ( i = iCD; i < itmp; i++ )
		{
			if ( lDisV[i] )
			{
				k--;
				for ( j = 0; j < N; j++ )
				{
					b[i][j] = b[k][j];
				}
			}
			if ( i > k )	break;
		}

		itmp = k;

		if (itmp == iCD)
		{
			Message( "What a strange matrix!" );
			return;
		}

		for ( i = iCD; i < itmp; i++ )
		{
			Hx(SA, HU, h[i], b[i]);
			iCount++;
		}

		for ( i = 0; i < iCD; i++ )
		{
			for ( k = iCD; k < itmp; k++ )
			{
				Ttmp[i*iDim+k] = 0.0;
				for ( j = 0; j < N; j++ )
					Ttmp[i*iDim+k] += h[i][j] * b[k][j];
				Ttmp[k*iDim+i] = Ttmp[i*iDim+k];
			}
		}

		for ( i = iCD; i < itmp; i++ )
		{
			for ( k = i; k < itmp; k++ )
			{
				Ttmp[i*iDim+k] = 0.0;
				for ( j = 0; j < N; j++ )
					Ttmp[i*iDim+k] += h[i][j] * b[k][j];
				Ttmp[k*iDim+i] = Ttmp[i*iDim+k];
			}
		}

		iCD = itmp;
	}

	for (i = 0; i < iEig; i++)
	{
		En[i] = Mu[i];
		for (j = 0; j < N; j++)
		{
			v[i * N + j] = 0.0;
			for (k = 0; k < iCD; k++)
				v[i * N + j] += b[k][j] * a[k][i];
		}
	}

	return;
}
