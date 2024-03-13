//
// File InfDMRGKondoQN.cpp
//
// Implementation of class CInfDMRGKondoQN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Kondo model
//
#include <iostream>
#include <vector>
#include <conio.h>
#include <time.h>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "CSparse.hpp"
#include "InfDMRGKondoQN.h"
#include "KondoOp.h"

//#include "Debug.hpp"

#define QUANTUM_NUMBERS	2

int N_L = 4;
int N_A = 4;
int N_LA = 16;
int N_d = 2;
int N_LA_d = 32;


Matrix<double> L2a(N_LA_d, N_LA_d, 0.0), A2a(N_LA_d, N_LA_d, 0.0);
Matrix<double> L3a(N_LA_d, N_LA_d, 0.0), A3a(N_LA_d, N_LA_d, 0.0);
Matrix<double> L2b(N_LA_d, N_LA_d, 0.0), A2b(N_LA_d, N_LA_d, 0.0);
Matrix<double> L3b(N_LA_d, N_LA_d, 0.0), A3b(N_LA_d, N_LA_d, 0.0);

Matrix<double> H2(N_LA_d, N_LA_d, 0.0), H3(N_LA_d, N_LA_d, 0.0);

Matrix<int>		QN_LA_2(2, N_LA_d), QN_LA_3(2, N_LA_d);
Vector<double>	QN_d_2 (N_LA_d),    QN_d_3 (N_LA_d);

vector<long, allocator<long> > QN_index_H;

vector<long, allocator<long> > vec_1;
vector<long, allocator<long> > vec_2;
vector<long, allocator<long> > vec_3;
vector<long, allocator<long> > vec_4;


CInfDMRGKondoQN::CInfDMRGKondoQN(int m, int BasisSize, double t, double alpha,
				double UL, double UA, double J, int Krylov)
{
	switch (BasisSize)
	{
	case 4:
		N_L = 2;
		N_A = 2;
		N_LA = 4;
		N_d = 1;
		N_LA_d = 4;
		break;
	case 18:
		N_L = 3;
		N_A = 3;
		N_LA = 9;
		N_d = 2;
		N_LA_d = 18;
		break;
	case 32:
		N_L = 4;
		N_A = 4;
		N_LA = 16;
		N_d = 2;
		N_LA_d = 32;
		break;
	default:
		Error("Wrong 'BasisSize' in 'CInfDMRGKondoQN'");
		break;
	}

	if (m < N_LA_d)	Error("Non-equivalent sizes in 'CInfDMRGKondoQN'");

	m_m = m;
	m_Size = m_m * N_LA_d * N_LA_d * m_m;
	m_EPS_DMRG = 1.0E-005;
	m_EPS_Diag = 1.0E-008;
	m_Iterations = 100;
	m_Krylov = Krylov;

	m_t = t;
	m_alpha = alpha;
	m_UL = UL;
	m_UA = UA;
	m_J = J;

	A1a = new Matrix<double>(m_m, m_m, 0.0);
	L4a = new Matrix<double>(m_m, m_m, 0.0);
	A1b = new Matrix<double>(m_m, m_m, 0.0);
	L4b = new Matrix<double>(m_m, m_m, 0.0);

	H1 = new Matrix<double>(m_m, m_m, 0.0);
	H4 = new Matrix<double>(m_m, m_m, 0.0);

	QN_LA_1 = new Matrix<int>(2, m_m * N_LA_d, 0);
	QN_d_1  = new Vector<double>(m_m * N_LA_d, 0.0);
	QN_LA_4 = new Matrix<int>(2, m_m * N_LA_d, 0);
	QN_d_4  = new Vector<double>(m_m * N_LA_d, 0.0);

	SetupStandardBlocks();
};

CInfDMRGKondoQN::~CInfDMRGKondoQN()
{
	delete QN_d_4;
	delete QN_LA_4;
	delete QN_d_1;
	delete QN_LA_1;

	delete H4;
	delete H1;

	delete L4b;
	delete A1b;
	delete L4a;
	delete A1a;
}

void CInfDMRGKondoQN::SetupStandardBlocks(void)
{
	for (int di = 0; di < N_d; di++)
	for (int ai = 0; ai < N_A; ai++)
	for (int li = 0; li < N_L; li++)
	{
		int L = di * N_LA + ai * N_L + li;

		double value = 0.0;

		value += m_alpha * (nLna[li] + nLnb[li]);
		value += m_UL * nLna[li] * nLnb[li];
		value += m_UA * nAna[ai] * nAnb[ai];
		value += m_J * SLn_z[li][li] * Sdn_z[di][di];

		(*H1)[L][L] = value;
		  H2 [L][L] = value;
		  H3 [L][L] = value;
		(*H4)[L][L] = value;

#if QUANTUM_NUMBERS == 3
		(*QN_LA_1)[0][L] = nAna[ai] + nLna[li];
		  QN_LA_2 [0][L] = nAna[ai] + nLna[li];
		  QN_LA_3 [0][L] = nAna[ai] + nLna[li];
		(*QN_LA_4)[0][L] = nAna[ai] + nLna[li];

		(*QN_LA_1)[1][L] = nAnb[ai] + nLnb[li];
		  QN_LA_2 [1][L] = nAnb[ai] + nLnb[li];
		  QN_LA_3 [1][L] = nAnb[ai] + nLnb[li];
		(*QN_LA_4)[1][L] = nAnb[ai] + nLnb[li];

		(*QN_d_1)[L] = 0.5 * (ndna[di] - ndnb[di]);
		  QN_d_2 [L] = 0.5 * (ndna[di] - ndnb[di]);
		  QN_d_3 [L] = 0.5 * (ndna[di] - ndnb[di]);
		(*QN_d_4)[L] = 0.5 * (ndna[di] - ndnb[di]);
#elif QUANTUM_NUMBERS == 2
		int n_up   = nLna[li] + nAna[ai] + ndna[di];
		int n_down = nLnb[li] + nAnb[ai] + ndnb[di];
		int n_LA   = nLna[li] + nAna[ai] + nLnb[li] + nAnb[ai];

		(*QN_LA_1)[0][L] = n_LA;
		  QN_LA_2 [0][L] = n_LA;
		  QN_LA_3 [0][L] = n_LA;
		(*QN_LA_4)[0][L] = n_LA;

		(*QN_LA_1)[1][L] = 0;
		  QN_LA_2 [1][L] = 0;
		  QN_LA_3 [1][L] = 0;
		(*QN_LA_4)[1][L] = 0;

		(*QN_d_1)[L] = 0.5 * (n_up - n_down);
		  QN_d_2 [L] = 0.5 * (n_up - n_down);
		  QN_d_3 [L] = 0.5 * (n_up - n_down);
		(*QN_d_4)[L] = 0.5 * (n_up - n_down);
#endif

		for (int dj = 0; dj < N_d; dj++)
		for (int aj = 0; aj < N_A; aj++)
		for (int lj = 0; lj < N_L; lj++)
		{
			int R = dj * N_LA + aj * N_L + lj;

			if (li == lj && di == dj)
			{
				A3a[L][R] = A2a[L][R] = (*A1a)[L][R] = Ana[ai][aj];
				A3b[L][R] = A2b[L][R] = (*A1b)[L][R] = Anb[ai][aj];
			}

			if (ai == aj && di == dj)
			{
				L3a[L][R] = L2a[L][R] = (*L4a)[L][R] = Lna[li][lj];
				L3b[L][R] = L2b[L][R] = (*L4b)[L][R] = Lnb[li][lj];
			}

			if (L > R || L == R)	continue;

			if (ai == aj)
			{
				value = 0.5 * m_J * (SLn_up[li][lj] * Sdn_up[dj][di] + SLn_up[lj][li] * Sdn_up[di][dj]);
				
				(*H1)[L][R] += value;
				  H2 [L][R] += value;
				  H3 [L][R] += value;
				(*H4)[L][R] += value;
			}

			if (di == dj)
			{
				value = 0.0;

				value += m_t * (Ana[lj][li] * Lna[ai][aj] + Lna[aj][ai] * Ana[li][lj]);
				value += m_t * (Anb[lj][li] * Lnb[ai][aj] + Lnb[aj][ai] * Anb[li][lj]);

				(*H1)[L][R] += value;
				  H2 [L][R] += value;
				  H3 [L][R] += value;
				(*H4)[L][R] += value;
			}

			(*H1)[R][L] = (*H1)[L][R];
			  H2 [R][L] =   H2 [L][R];
			  H3 [R][L] =   H3 [L][R];
			(*H4)[R][L] = (*H4)[L][R];
		}
	}
}

int CInfDMRGKondoQN::GetSubSize(int* mi, int qn_LA, int qn_d, double qn_S)
{
	int _m234 = mi[1] * mi[2] * mi[3];
	int _m34 = mi[2] * mi[3];
	int _m4 = mi[3];

	int Counter = 0;

	QN_index_H.resize(0);

	for (int k1 = 0; k1 < mi[0]; k1++)
	for (int k2 = 0; k2 < mi[1]; k2++)
	for (int k3 = 0; k3 < mi[2]; k3++)
	for (int k4 = 0; k4 < mi[3]; k4++)
	{
		long L = k1 * _m234 + k2 * _m34 + k3 * _m4 + k4;

#if QUANTUM_NUMBERS == 3
		int QN_LA_up   = (*QN_LA_1)[0][k1] + QN_LA_2[0][k2] + QN_LA_3[0][k3] + (*QN_LA_4)[0][k4];
		int QN_LA_down = (*QN_LA_1)[1][k1] + QN_LA_2[1][k2] + QN_LA_3[1][k3] + (*QN_LA_4)[1][k4];
		double QN_d = (*QN_d_1)[k1] + QN_d_2[k2] + QN_d_3[k3] + (*QN_d_4)[k4];
		double QN_S = 0.5 * (QN_LA_up - QN_LA_down) + QN_d;

		if ((QN_LA_up + QN_LA_down) == qn_LA && QN_S == qn_S)
		{
			QN_index_H.push_back(L);
			Counter++;
		}
#elif QUANTUM_NUMBERS == 2
		int QN_LA = (*QN_LA_1)[0][k1] + QN_LA_2[0][k2] + QN_LA_3[0][k3] + (*QN_LA_4)[0][k4];
		double QN_S = (*QN_d_1)[k1] + QN_d_2[k2] + QN_d_3[k3] + (*QN_d_4)[k4];

		if (QN_LA == qn_LA && QN_S == qn_S)
		{
			QN_index_H.push_back(L);
			Counter++;
		}
#endif

	}

	return Counter;
}

int CInfDMRGKondoQN::IsExit(int Param, double Convergence, int Iteration)
{
	if (Param == 0)
	{
		return (fabs(Convergence) > m_EPS_DMRG);
	}
	else if (Param == 1)
	{
		return (Iteration < m_Iterations);
	}
	else
	{
		Error("Wrong value of 'Param'");
		return 0;
	}
}

double CInfDMRGKondoQN::Iterations(int& L, Vector<double>& TargetState, int qn_LA, int qn_d, double qn_S, int Param)
{
	int m1 = N_LA_d, m2 = N_LA_d;
	int m1m2 = m_m * m2;
	int mi[4] = { m1, m2, m2, m1 };
	int SitesNumber = 4;

	cout << "m1 = " << m1 << endl;
	cout << "m1m2 = " << m1m2 << endl;

	double	GSE,
			GSE_per_site = 0.0,
			Convergence = 1.0,
			MaxAccuracy = 0.0;

	long	MaxSize = TargetState.getSize();

	Vector<double> DM12(m1m2 * m1m2, 0.0);
	Vector<double> w1(m1m2, 0.0);
	Vector<double> O1(m1m2 * m1m2, 0.0);

	cout << "Memory allocated successfully..." << endl;

	Matrix<int>		QN_LA(2, m1m2, 0);
	Vector<double>	QN_d(m1m2, 0.0);
	Vector<long>	Multi_Index(m1m2);

	CSparseSymmetric<double, long> SubHSB;
	CSparseSymmetric<double, long> H12;
	CSparseUnitary<double, long> __A1a, __L4a;
	CSparseUnitary<double, long> __A1b, __L4b;

	struct tm *newtime;
	time_t start, finish;

	time(&start);
	newtime = localtime(&start);

	FILE* stream = fopen("DAT\\InfDMRGKondoQN.dat", "w+");
		fprintf(stream, "DMRG RESULTS:\n");
		fprintf(stream, "\nSTART: %s", asctime(newtime));
	fclose(stream);

	int Iteration = 0;

	while (IsExit(Param, Convergence, Iteration))
	{
		++Iteration;
		L = SitesNumber;
		m1m2 = m1 * m2;
		m_Size = m1 * m2 * m2 * m1;

		int		_QN_LA = SitesNumber * qn_LA;
		int		_QN_d = SitesNumber * qn_d;
		double	_QN_S;
		
		if (qn_S == 0 || qn_S == 4)
		{
			_QN_S = SitesNumber * qn_S / 4;
		}
		else
		{
			_QN_S = SitesNumber - (4 - qn_S);
		}

		int SubSize = GetSubSize(mi, _QN_LA, _QN_d, _QN_S);

		if (SubSize > MaxSize)
		{
			MaxSize = SubSize;
			TargetState.Realloc(SubSize);
		}
		
		SubHSB.reserve(SubSize, SubSize, SubSize);
		SubHSB.resize(SubSize, 0, 0);

		cout << "\nIteration = " << Iteration << endl;
		cout << "QN_LA = " << _QN_LA << "; QN_d = " << _QN_d << "; QN_S = " << _QN_S << endl;
		cout << "mi : " << mi[0] << ", " << mi[1] << ", " << mi[2] << ", " << mi[3] << endl;
		cout << SitesNumber << " sites" << endl;
		cout << "Size = " << m_Size;
		cout << ";\t SubSize = " << SubSize << endl;

		FormSuperBlock(mi, SubHSB);
		
		cout << "Sub : NEL = " << SubHSB.ir.size();
		cout << ";\t Sparsity = " << SubHSB.GetSparsity() << " %" << endl;

		GSE = Diagonalize_HSB(SubSize, SubHSB, TargetState);

		Convergence = GSE_per_site - GSE / SitesNumber;
		GSE_per_site = GSE / SitesNumber;

		cout << "Sub : GSE = " << GSE;
		cout << "; GSE / L = " << GSE_per_site;
		cout << "; Convergence = " << Convergence << endl;

		stream = fopen("DAT\\InfDMRGKondoQN.dat", "a+");
			fprintf(stream, "\nE(L=%d,\tm=%d)", SitesNumber, mi[0]);
			fprintf(stream, "\t= %lf;", GSE);
			fprintf(stream, "\tE / L = %lf;", GSE_per_site);
			fprintf(stream, "\tTruncation error = %6.5e;", MaxAccuracy);
			fprintf(stream, "\tIteration = %d;", Iteration);
			fprintf(stream, "\tConvergence = %6.5e;", Convergence);
		fclose(stream);

		int _m1 = (m1m2 > m_m) ? m_m : m1m2;

//////////////////////////////////////////////////////////////////////////////////

		FormDensityMatrix12(mi, DM12, TargetState, QN_LA, QN_d);
		Arrange_ro(m1m2, QN_LA, QN_d, *QN_LA_1, *QN_d_1, DM12, Multi_Index);

		double Accuracy_1 = Diagonalize_ro(m1m2, DM12, w1, O1, QN_LA, QN_d, Multi_Index);
		Form12BlockSystem(mi, H12, __A1a, __A1b, Multi_Index);
		Arrange_QN(_m1, m1m2, O1, QN_LA, QN_d, *QN_LA_1, *QN_d_1);

		Truncation(_m1, m1m2, DM12, H12,   O1, *H1);
		Truncation(_m1, m1m2, DM12, __A1a, O1, *A1a);
		Truncation(_m1, m1m2, DM12, __A1b, O1, *A1b);

//////////////////////////////////////////////////////////////////////////////////

		FormDensityMatrix34(mi, DM12, TargetState, QN_LA, QN_d);
		Arrange_ro(m1m2, QN_LA, QN_d, *QN_LA_4, *QN_d_4, DM12, Multi_Index);

		double Accuracy_4 = Diagonalize_ro(m1m2, DM12, w1, O1, QN_LA, QN_d, Multi_Index);
		Form34BlockSystem(mi, H12, __L4a, __L4b, Multi_Index);
		Arrange_QN(_m1, m1m2, O1, QN_LA, QN_d, *QN_LA_4, *QN_d_4);

		Truncation(_m1, m1m2, DM12, H12,   O1, *H4);
		Truncation(_m1, m1m2, DM12, __L4a, O1, *L4a);
		Truncation(_m1, m1m2, DM12, __L4b, O1, *L4b);

//////////////////////////////////////////////////////////////////////////////////

		m1 = (m1m2 > m_m) ? m_m : m1m2;
		mi[0] = mi[3] = m1;

		MaxAccuracy = (fabs(Accuracy_1) > fabs(Accuracy_4)) ? Accuracy_1 : Accuracy_4;

		SitesNumber += 2;
	}

	time(&finish);
	newtime = localtime(&finish);
	double elapsed_time = difftime(finish, start);
	
	int hours = (int)(elapsed_time / 3600.0);
	int minutes = (int)((elapsed_time - hours * 3600) / 60.0);
	int seconds = (int)(elapsed_time - hours * 3600 - minutes * 60);

	stream = fopen("DAT\\InfDMRGKondoQN.dat", "a+");
		fprintf(stream, "\n\nFINISH: %s", asctime(newtime));
		fprintf(stream, "\nELAPSED TIME: %d hours %d minutes %d seconds",
				hours, minutes, seconds);
	fclose(stream);

	return GSE;
}

void CInfDMRGKondoQN::FormSuperBlock(int* mi, CSparseSymmetric<double, long>& SubHSB)
{
	int _m234 = mi[1] * mi[2] * mi[3];
	int _m34 = mi[2] * mi[3];
	int _m4 = mi[3];

	for (int Index_L = 0; Index_L < QN_index_H.size(); Index_L++)
	{
		long L = QN_index_H[Index_L];

		int i1 = L / _m234;
		int i2_i3_i4 = L % _m234;
		int i2 = i2_i3_i4 / _m34;
		int i3_i4 = i2_i3_i4 % _m34;
		int i3 = i3_i4 / _m4;
		int i4 = i3_i4 % _m4;
		
		for (int Index_R = Index_L; Index_R < QN_index_H.size(); Index_R++)
		{
			long R = QN_index_H[Index_R];

			int j1 = R / _m234;
			int j2_j3_j4 = R % _m234;
			int j2 = j2_j3_j4 / _m34;
			int j3_j4 = j2_j3_j4 % _m34;
			int j3 = j3_j4 / _m4;
			int j4 = j3_j4 % _m4;

			double Summa_B = 0.0;

			if (i2 == j2 && i3 == j3 && i4 == j4)
			{
				Summa_B += (*H1)[i1][j1];
			}

			if (i1 == j1 && i3 == j3 && i4 == j4)
			{
				Summa_B += H2[i2][j2];
			}

			if (i1 == j1 && i2 == j2 && i4 == j4)
			{
				Summa_B += H3[i3][j3];
			}

			if (i1 == j1 && i2 == j2 && i3 == j3)
			{
				Summa_B += (*H4)[i4][j4];
			}

			double Summa_t = 0.0;

			if (i3 == j3 && i4 == j4)
			{
				Summa_t += (*A1a)[j1][i1] * L2a[i2][j2] + L2a[j2][i2] * (*A1a)[i1][j1];
				Summa_t += (*A1b)[j1][i1] * L2b[i2][j2] + L2b[j2][i2] * (*A1b)[i1][j1];
			}

			if (i1 == j1 && i4 == j4)
			{
				Summa_t += A2a[j2][i2] * L3a[i3][j3] + L3a[j3][i3] * A2a[i2][j2];
				Summa_t += A2b[j2][i2] * L3b[i3][j3] + L3b[j3][i3] * A2b[i2][j2];
			}

			if (i1 == j1 && i2 == j2)
			{
				Summa_t += A3a[j3][i3] * (*L4a)[i4][j4] + (*L4a)[j4][i4] * A3a[i3][j3];
				Summa_t += A3b[j3][i3] * (*L4b)[i4][j4] + (*L4b)[j4][i4] * A3b[i3][j3];
			}

			double Summa = Summa_t * m_t + Summa_B;
			
			SubHSB.push_back(Summa, Index_L, Index_R);
		}
	}
}

double CInfDMRGKondoQN::Diagonalize_HSB(int Size, CSparseSymmetric<double, long>& HSB, Vector<double>& TargetState)
{
	double E = 0.0;

	if (TargetState.getSize() < Size)
	{
		Error("Size of TargetState is too small in 'Diagonalize_HSB'");
	}

	if (Size < 300)
	{
		Matrix<double> A2(Size, Size, 0.0);
		Vector<double> d(Size+1, 0.0);
		Vector<double> e(Size+1, 0.0);

		for (int ii = 0; ii < Size; ii++)
		{
			A2[ii][ii] = HSB.sa[ii];
		}

		for (int ii = 0; ii < HSB.ir.size(); ii++)
		{
			int i = HSB.ir[ii];
			int j = HSB.ic[ii];

			A2[j][i] = A2[i][j] = HSB.sa[ii + Size];
		}

		int SMALL = 1;

		tred2_c(A2.getMat(), Size, d.getVec(), e.getVec());
		tqli_c(d.getVec(), e.getVec(), Size, A2.getMat());
		sort_c(Size, SMALL, d.getVec(), A2.getMat());

		E = d[0];

		for (int ii = 0; ii < Size; ii++)
		{
			TargetState[ii] = A2[ii][0];
		}
	}
	else
	{
		SparseMatrix<double> SA(Size, HSB.ir.size(), HSB.sa.data(),
			HSB.ir.data(), HSB.ic.data(), SYMMETRIC | FULLDIAG);

		int		NEVEC = 1;
		int		nDim = m_Krylov * NEVEC;

		Vector<double>	EVAL(NEVEC);

		Davidson(NEVEC, SA, EVAL, TargetState, nDim, m_EPS_Diag);

		E = EVAL[0];
	}

	return E;
}

void CInfDMRGKondoQN::FormDensityMatrix12(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Matrix<int>& QN_LA, Vector<double>& QN_d)
{
	int _m234 = mi[1] * mi[2] * mi[3];
	int _m12 = mi[0] * mi[1];
	int _m34 = mi[2] * mi[3];
	int _m4 = mi[3];

	int SubSize = QN_index_H.size();

	if (_m12 * _m12 > DM.getSize())
	{
		Error("Non-equivalent sizes in 'FormDensityMatrix12'");
	}

	vec_1.reserve(SubSize);
	vec_2.reserve(SubSize);
	vec_3.reserve(SubSize);
	vec_4.reserve(SubSize);

	vec_1.resize(SubSize);
	vec_2.resize(SubSize);
	vec_3.resize(SubSize);
	vec_4.resize(SubSize);

	for (int index = 0; index < SubSize; index++)
	{
		long L = QN_index_H[index];

		int i1 = L / _m234;
		int i2_i3_i4 = L % _m234;
		int i2 = i2_i3_i4 / _m34;
		int i3_i4 = i2_i3_i4 % _m34;
		int i3 = i3_i4 / _m4;
		int i4 = i3_i4 % _m4;

		vec_1[index] = i1;
		vec_2[index] = i2;
		vec_3[index] = i3;
		vec_4[index] = i4;
	}

	for (int i = 0; i < _m12; i++)
	{
		DM[i * _m12 + i] = 0.0;

		int i1 = i / mi[1];
		int i2 = i % mi[1];

		QN_LA[0][i] = (*QN_LA_1)[0][i1] + QN_LA_2[0][i2];
		QN_LA[1][i] = (*QN_LA_1)[1][i1] + QN_LA_2[1][i2];

		QN_d[i] = (*QN_d_1)[i1] + QN_d_2[i2];

		for (int j = i + 1; j < _m12; j++)
		{
			DM[j * _m12 + i] = DM[i * _m12 + j] = 0.0;
		}
	}

	for (int i_Index = 0; i_Index < SubSize;)
	{
		int i1 = vec_1[i_Index];
		int i2 = vec_2[i_Index];

		for (int j_Index = i_Index; j_Index < SubSize;)
		{
			int j1 = vec_1[j_Index];
			int j2 = vec_2[j_Index];

			long LL = i1 * mi[1] + i2;
			long RR = j1 * mi[1] + j2;

			double Summa = 0.0;
			
			int i_Summator = i_Index;
			int j_Summator = j_Index;

			for (int i3 = vec_3[i_Summator]; i1 == vec_1[i_Summator] && i2 == vec_2[i_Summator];)
			{
				int i4 = vec_4[i_Summator];

				for (int j3 = vec_3[j_Summator]; j1 == vec_1[j_Summator] && j2 == vec_2[j_Summator];)
				{
					int j4 = vec_4[j_Summator];

					if (j3 < i3)
					{
						j3 = vec_3[++j_Summator];
					}
					else if (j3 == i3)
					{
						if (j4 < i4)
						{
							j3 = vec_3[++j_Summator];
						}
						else if (j4 == i4)
						{
							Summa += TargetState[i_Summator] * TargetState[j_Summator];
							j_Summator++;
							break;
						}
						else if (j4 > i4)
						{
							break;
						}
					}
					else if (j3 > i3)
					{
						break;
					}
				}
				
				i3 = vec_3[++i_Summator];
			}

			DM[RR * _m12 + LL] = DM[LL * _m12 + RR] = Summa;

			while (++j_Index < SubSize)
			{
				if (j1 != vec_1[j_Index] || j2 != vec_2[j_Index])	break;
			}
		}

		while (++i_Index < SubSize)
		{
			if (i1 != vec_1[i_Index] || i2 != vec_2[i_Index])	break;
		}
	}
}

void CInfDMRGKondoQN::FormDensityMatrix34(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Matrix<int>& QN_LA, Vector<double>& QN_d)
{
	int _m34 = mi[2] * mi[3];

	int SubSize = QN_index_H.size();

	for (int i = 0; i < SubSize - 1; i++)
	{
		int k = i;

		int k3 = vec_3[i];
		int k4 = vec_4[i];
		int k1 = vec_1[i];
		int k2 = vec_2[i];

		for (int j = i + 1; j < SubSize; j++)
		{
			if (vec_3[j] < k3 || (vec_3[j] == k3 && vec_4[j] < k4) ||
				(vec_3[j] == k3 && vec_4[j] == k4 && vec_1[j] < k1) ||
				(vec_3[j] == k3 && vec_4[j] == k4 && vec_1[j] == k1 && vec_2[j] < k2))
			{
				k = j;

				k3 = vec_3[j];
				k4 = vec_4[j];
				k1 = vec_1[j];
				k2 = vec_2[j];
			}
		}

		if (k != i)
		{
			vec_3[k] = vec_3[i];
			vec_3[i] = k3;

			vec_4[k] = vec_4[i];
			vec_4[i] = k4;

			vec_1[k] = vec_1[i];
			vec_1[i] = k1;

			vec_2[k] = vec_2[i];
			vec_2[i] = k2;

			int L = QN_index_H[i];
			double Value = TargetState[i];

			QN_index_H[i] = QN_index_H[k];
			TargetState[i] = TargetState[k];

			QN_index_H[k] = L;
			TargetState[k] = Value;
		}
	}

	for (int i = 0; i < _m34; i++)
	{
		DM[i * _m34 + i] = 0.0;

		int i3 = i / mi[3];
		int i4 = i % mi[3];

		QN_LA[0][i] = QN_LA_3[0][i3] + (*QN_LA_4)[0][i4];
		QN_LA[1][i] = QN_LA_3[1][i3] + (*QN_LA_4)[1][i4];
		
		QN_d[i] = QN_d_3[i3] + (*QN_d_4)[i4];

		for (int j = i + 1; j < _m34; j++)
		{
			DM[j * _m34 + i] = DM[i * _m34 + j] = 0.0;
		}
	}

	for (int i_Index = 0; i_Index < SubSize;)
	{
		int i3 = vec_3[i_Index];
		int i4 = vec_4[i_Index];

		for (int j_Index = i_Index; j_Index < SubSize;)
		{
			int j3 = vec_3[j_Index];
			int j4 = vec_4[j_Index];

			long LL = i3 * mi[3] + i4;
			long RR = j3 * mi[3] + j4;

			double Summa = 0.0;
			
			int i_Summator = i_Index;
			int j_Summator = j_Index;

			for (int i1 = vec_1[i_Summator]; i3 == vec_3[i_Summator] && i4 == vec_4[i_Summator];)
			{
				int i2 = vec_2[i_Summator];

				for (int j1 = vec_1[j_Summator]; j3 == vec_3[j_Summator] && j4 == vec_4[j_Summator];)
				{
					int j2 = vec_2[j_Summator];

					if (j1 < i1)
					{
						j1 = vec_1[++j_Summator];
					}
					else if (j1 == i1)
					{
						if (j2 < i2)
						{
							j1 = vec_1[++j_Summator];
						}
						else if (j2 == i2)
						{
							Summa += TargetState[i_Summator] * TargetState[j_Summator];
							j_Summator++;
							break;
						}
						else if (j2 > i2)
						{
							break;
						}
					}
					else if (j1 > i1)
					{
						break;
					}
				}
				
				i1 = vec_1[++i_Summator];
			}

			DM[RR * _m34 + LL] = DM[LL * _m34 + RR] = Summa;

			while (++j_Index < SubSize)
			{
				if (j3 != vec_3[j_Index] || j4 != vec_4[j_Index])	break;
			}
		}

		while (++i_Index < SubSize)
		{
			if (i3 != vec_3[i_Index] || i4 != vec_4[i_Index])	break;
		}
	}
}

double CInfDMRGKondoQN::Diagonalize_ro(int Size, Vector<double>& DM, Vector<double>& w, Vector<double>& O,
					Matrix<int>& QN_LA, Vector<double>& QN_d, Vector<long>& Multi_Index)
{
	int __m = (Size > m_m) ? m_m : Size;

	if (DM.getSize() < Size * Size || w.getSize() < Size ||
		O.getSize() < Size * Size)
	{
		Error("Wrong dimensions in 'Diagonalize_ro'");
	}

	int MaxSubSize = 1;
	int SubSize = 0;
	int SubNumber = 1;

	Vector<int> SubSpaces(Size, 0);

	for (int i = 0; i < Size; i++)
	{
		SubSize++;

		if (i == Size - 1)
		{
			SubSpaces[SubNumber - 1] = SubSize;
			break;
		}

		if (QN_LA[0][i] != QN_LA[0][i + 1] || QN_d[i] != QN_d[i + 1] ||
			QN_LA[1][i] != QN_LA[1][i + 1])
		{
			SubSpaces[SubNumber - 1] = SubSize;
			SubNumber++;
			if (SubSize > MaxSubSize)	MaxSubSize = SubSize;
			SubSize = 0;
		}
	}

	Vector<double> Sub_DM(MaxSubSize * MaxSubSize, 0.0);
	Vector<double> Sub_w(MaxSubSize, 0.0);
	Vector<double> Sub_O(MaxSubSize * MaxSubSize, 0.0);

	for (int i = 0; i < Size; i++)
	{
		O[i * Size + i] = w[i] = 0.0;

		for (int j = i + 1; j < Size; j++)
		{
			O[i * Size + j] = O[j * Size + i] = 0.0;
		}
	}

	int Shift = 0;

	for (int i = 0; i < SubNumber; i++)
	{
		SubSize = SubSpaces[i];

		for (int ii = 0; ii < SubSize; ii++)
		{
			Sub_w[ii] = 0.0;

			for (int jj = 0; jj < SubSize; jj++)
			{
				Sub_O[ii * SubSize + jj] = 0.0;
				Sub_DM[ii * SubSize + jj] = DM[(Shift + ii) * Size + Shift + jj];
			}
		}

		if (SubSize > 1)
		{
			Matrix<double> A2(SubSize, SubSize, 0.0);
        	        Vector<double> d(SubSize+1, 0.0);
                	Vector<double> e(SubSize+1, 0.0);

	                for (int i = 0; i < SubSize; i++)
        	        {
                	   for (int j = 0; j < SubSize; j++)
	                   {
			      A2[i][j] = Sub_DM[i * SubSize + j];
	                   }
	                }

			tred2_c(A2.getMat(), SubSize, d.getVec(), e.getVec());
			tqli_c(d.getVec(), e.getVec(), SubSize, A2.getMat());

	                for (int i = 0; i < SubSize; i++)
	                {
			   Sub_w[i] = d[i];

	                   for (int j = 0; j < SubSize; j++)
	                   {
			      Sub_DM[i * SubSize + j] = A2[j][i];
	                   }
	                }

                        for (int jj = 0; jj < SubSize * SubSize; jj++)
                        {
                                Sub_O[jj] = Sub_DM[jj];
                        }
		}
		else if (SubSize == 1)
		{
			Sub_w[0] = Sub_DM[0];
			Sub_O[0] = 1.0;
		}
		else
		{
			Error("Something is wrong with SubSpaces in 'Diagonalize_ro'");
		}

		for (int ii = 0; ii < SubSize; ii++)
		{
			w[Shift + ii] = Sub_w[ii];

			for (int jj = 0; jj < SubSize; jj++)
			{
				O[(Shift + ii) * Size + Shift + jj] = Sub_O[ii * SubSize + jj];
			}
		}

		Shift += SubSpaces[i];
	}

	int SMALL = 0;

        Vector<int> index(Size, 0);

        for (int k1 = 0; k1 < Size; k1++)
        {
           index[k1] = k1;
        }

	sort_qn_3_c(Size, SMALL, w.getVec(), O.getVec(), QN_LA.getMat(), QN_d.getVec(), Multi_Index.getVec(), index.getVec());

	double P = 0.0;
	double Spur_ro = 0.0;
	double Spur_O = 0.0;

	for (int i = 0; i < __m; i++)
	{
		P += w[i];
	}

	for (int i = 0; i < Size; i++)
	{
		Spur_ro += w[i];

		for (int k = 0; k < Size; k++)
		{
			Spur_O += O[i * Size + k] * O[i * Size + k];
		}
	}

	double Accuracy = 1.0 - P;

	cout << "Spur(O) = " << Spur_O;
	cout << ";\t Spur(ro) = " << Spur_ro << ";\t Truncation error = " << Accuracy << endl;
												 
	if (fabs(Spur_O - Size) > 1.0E-008 || fabs(Spur_ro - 1) > 1.0E-008)
	{
		cout << "Spur(O) - Size = " << (Spur_O - Size);
		cout << ";\t Spur(ro) - 1 = " << (Spur_ro - 1) << endl;
		if (_getche() == 3)	exit(0);
	}

	return Accuracy;
}

void CInfDMRGKondoQN::Form12BlockSystem(int* mi, CSparseSymmetric<double, long>& H12,
					CSparseUnitary<double, long>& __A1a,
					CSparseUnitary<double, long>& __A1b,
					Vector<long>& Multi_Index)
{
	int m1m2 = mi[0] * mi[1];

	H12.reserve(m1m2, m1m2, m1m2);
	H12.resize(m1m2, 0, 0);

	__A1a.resize(m1m2 + 1, m1m2 + 1);
	__A1a.ija[0] = m1m2 + 1;

	__A1b.resize(m1m2 + 1, m1m2 + 1);
	__A1b.ija[0] = m1m2 + 1;

	for (int i = 0; i < m1m2; i++)
	{
		long L = Multi_Index[i];
		int i1 = L / mi[1];
		int i2 = L % mi[1];

		__A1a.sa[i] = 0.0;
		__A1b.sa[i] = 0.0;

		for (int j = 0; j < m1m2; j++)
		{
			long R = Multi_Index[j];
			int j1 = R / mi[1];
			int j2 = R % mi[1];
		
			if (i1 == j1)
			{
				__A1a.push_back(A2a[i2][j2], i, j);
				__A1b.push_back(A2b[i2][j2], i, j);
			}
			
			if (i > j)	continue;

			double Summa_B = 0.0;

			if (i2 == j2)
			{
				Summa_B += (*H1)[i1][j1];
			}

			if (i1 == j1)
			{
				Summa_B += H2[i2][j2];
			}

			double Summa_t = 0.0;

			Summa_t += (*A1a)[j1][i1] * L2a[i2][j2] + L2a[j2][i2] * (*A1a)[i1][j1];
			Summa_t += (*A1b)[j1][i1] * L2b[i2][j2] + L2b[j2][i2] * (*A1b)[i1][j1];

			double Summa = Summa_t * m_t + Summa_B;

			H12.push_back(Summa, i, j);
		}

		__A1a.ija[i + 1] = __A1a.ija.size();
		__A1b.ija[i + 1] = __A1b.ija.size();
	}
}

void CInfDMRGKondoQN::Form34BlockSystem(int* mi, CSparseSymmetric<double, long>& H34,
					CSparseUnitary<double, long>& __L4a,
					CSparseUnitary<double, long>& __L4b,
					Vector<long>& Multi_Index)
{
	int m3m4 = mi[2] * mi[3];

	H34.reserve(m3m4, m3m4, m3m4);
	H34.resize(m3m4, 0, 0);

	__L4a.resize(m3m4 + 1, m3m4 + 1);
	__L4a.ija[0] = m3m4 + 1;

	__L4b.resize(m3m4 + 1, m3m4 + 1);
	__L4b.ija[0] = m3m4 + 1;

	for (int i = 0; i < m3m4; i++)
	{
		long L = Multi_Index[i];
		int i3 = L / mi[3];
		int i4 = L % mi[3];

		__L4a.sa[i] = 0.0;
		__L4b.sa[i] = 0.0;

		for (int j = 0; j < m3m4; j++)
		{
			long R = Multi_Index[j];
			int j3 = R / mi[3];
			int j4 = R % mi[3];

			if (i4 == j4)
			{
				__L4a.push_back(L3a[i3][j3], i, j);
				__L4b.push_back(L3b[i3][j3], i, j);
			}

			if (i > j)	continue;

			double Summa_B = 0.0;

			if (i4 == j4)
			{
				Summa_B += H3[i3][j3];
			}

			if (i3 == j3)
			{
				Summa_B += (*H4)[i4][j4];
			}

			double Summa_t = 0.0;

			Summa_t += A3a[j3][i3] * (*L4a)[i4][j4] + (*L4a)[j4][i4] * A3a[i3][j3];
			Summa_t += A3b[j3][i3] * (*L4b)[i4][j4] + (*L4b)[j4][i4] * A3b[i3][j3];

			double Summa = Summa_t * m_t + Summa_B;

			H34.push_back(Summa, i, j);
		}

		__L4a.ija[i + 1] = __L4a.ija.size();
		__L4b.ija[i + 1] = __L4b.ija.size();
	}
}

void CInfDMRGKondoQN::Arrange_ro(int m1m2, Matrix<int>& Source_LA, Vector<double>& Source_d,
				Matrix<int>& Destin_LA, Vector<double>& Destin_d,
				Vector<double>& DM, Vector<long>& Multi_Index)
{
	for (int i = 0; i < m1m2; i++)
	{
		Destin_LA[0][i] = Source_LA[0][i];
		Destin_LA[1][i] = Source_LA[1][i];
		Destin_d[i] = Source_d[i];

		Multi_Index[i] = i;
	}

	for (int i = 0; i < m1m2 - 1; i++)
	{
		int k = i;

		int qn_LA_up = Source_LA[0][i];
		int qn_LA_down = Source_LA[1][i];
		double qn_d = Source_d[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (Source_LA[0][j] > qn_LA_up ||
				(Source_LA[0][j] == qn_LA_up && Source_LA[1][j] > qn_LA_down) ||
				(Source_LA[0][j] == qn_LA_up && Source_LA[1][j] == qn_LA_down && Source_d[j] > qn_d))
			{
				k = j;
				
				qn_LA_up = Source_LA[0][j];
				qn_LA_down = Source_LA[1][j];
				qn_d = Source_d[j];
			}
		}

		if (k != i)
		{
			Source_LA[0][k] = Source_LA[0][i];
			Source_LA[0][i] = qn_LA_up;

			Source_LA[1][k] = Source_LA[1][i];
			Source_LA[1][i] = qn_LA_down;

			Source_d[k] = Source_d[i];
			Source_d[i] = qn_d;

			long L = Multi_Index[i];
			Multi_Index[i] = Multi_Index[k];
			Multi_Index[k] = L;

			for (int j = 0; j < m1m2; j++)
			{
				qn_d = DM[j * m1m2 + i];
				DM[j * m1m2 + i] = DM[j * m1m2 + k];
				DM[j * m1m2 + k] = qn_d;
			}
		}
	}

	for (int i = 0; i < m1m2 - 1; i++)
	{
		int k = i;

		int qn_LA_up = Destin_LA[0][i];
		int qn_LA_down = Destin_LA[1][i];
		double qn_d = Destin_d[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (Destin_LA[0][j] > qn_LA_up ||
				(Destin_LA[0][j] == qn_LA_up && Destin_LA[1][j] > qn_LA_down) ||
				(Destin_LA[0][j] == qn_LA_up && Destin_LA[1][j] == qn_LA_down && Destin_d[j] > qn_d))
			{
				k = j;
				
				qn_LA_up = Destin_LA[0][j];
				qn_LA_down = Destin_LA[1][j];
				qn_d = Destin_d[j];
			}
		}

		if (k != i)
		{
			Destin_LA[0][k] = Destin_LA[0][i];
			Destin_LA[0][i] = qn_LA_up;

			Destin_LA[1][k] = Destin_LA[1][i];
			Destin_LA[1][i] = qn_LA_down;

			Destin_d[k] = Destin_d[i];
			Destin_d[i] = qn_d;

			for (int j = 0; j < m1m2; j++)
			{
				qn_d = DM[i * m1m2 + j];
				DM[i * m1m2 + j] = DM[k * m1m2 + j];
				DM[k * m1m2 + j] = qn_d;
			}
		}
	}
}

void CInfDMRGKondoQN::Arrange_QN(int m, int m1m2, Vector<double>& O, Matrix<int>& Source_LA, Vector<double>& Source_d,
				Matrix<int>& Destin_LA, Vector<double>& Destin_d)
{
	for (int i = 0; i < m; i++)
	{
		Destin_LA[0][i] = Source_LA[0][i];
		Destin_LA[1][i] = Source_LA[1][i];
		Destin_d[i] = Source_d[i];
	}
}

void CInfDMRGKondoQN::Truncation(int m, int m1m2, Vector<double>& H12_O,
								 CSparse<double, long>& H12,
								 Vector<double>& O, Matrix<double>& H1)
{
	H12.Multiply_TR(m1m2, m, H12_O, O);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			H1[i][j] = 0.0;

			for (int k = 0; k < m1m2; k++)
			{
				H1[i][j] += O[i * m1m2 + k] * H12_O[k * m + j];
			}
		}
	}
}
