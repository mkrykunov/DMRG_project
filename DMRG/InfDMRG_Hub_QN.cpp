//
// File InfDMRG_Hub_QN.cpp
//
// Implementation of class CInfDMRG_Hub_QN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Hubbard model
//
#include <iostream>
#include <vector>
#include <conio.h>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "HubbardOp.h"
#include "CSparse.hpp"
#include "InfDMRG_Hub_QN.h"

#include "Debug.hpp"

const int N = 4;

Matrix<double> c2a(N, N), c3a(N, N);
Matrix<double> c2b(N, N), c3b(N, N);

Matrix<double> H2(N, N, 0.0), H3(N, N, 0.0);

Vector<int> QN_up_2(N), QN_up_3(N);
Vector<int> QN_down_2(N), QN_down_3(N);

vector<long, allocator<long> > QN_index_H;

vector<long, allocator<long> > vec_1;
vector<long, allocator<long> > vec_2;
vector<long, allocator<long> > vec_3;
vector<long, allocator<long> > vec_4;


CInfDMRG_Hub_QN::CInfDMRG_Hub_QN(int m, double t1, double t2, double U)
{
	m_m = m;
	m_t1 = t1;
	m_t2 = t2;
	m_Size = m_m * N * N * m_m;
	m_Accuracy = 0.0;
	m_EPS_DMRG = 1.0E-005;
	m_EPS_Diag = 1.0E-008;

	c1a = new Matrix<double>(m_m, m_m);
	c4a = new Matrix<double>(m_m, m_m);
	c1b = new Matrix<double>(m_m, m_m);
	c4b = new Matrix<double>(m_m, m_m);

	H1 = new Matrix<double>(m_m, m_m, 0.0);
	H4 = new Matrix<double>(m_m, m_m, 0.0);

	QN_up_1 = new Vector<int>(m_m * N, 0);
	QN_down_1 = new Vector<int>(m_m * N, 0);
	QN_up_4 = new Vector<int>(m_m * N, 0);
	QN_down_4 = new Vector<int>(m_m * N, 0);

	Initialization(U);
};

CInfDMRG_Hub_QN::~CInfDMRG_Hub_QN()
{
	delete c1a;
	delete c4a;
	delete c1b;
	delete c4b;

	delete H1;
	delete H4;

	delete QN_up_1;
	delete QN_down_1;
	delete QN_up_4;
	delete QN_down_4;
}

double CInfDMRG_Hub_QN::Iterations(int& L, Vector<double>& TargetState, int qn_up, int qn_down)
{
	int m1 = N, m2 = N;
	int m1m2 = m_m * m2;
	int mi[4] = { m1, m2, m2, m1 };
	int SitesNumber = 4;

	double GSE, GSE_prev = 0.0, Convergence = 1.0;
	long MaxSize = TargetState.getSize();

	Vector<double> DM12(m1m2 * m1m2, 0.0);
	Vector<double> w1(m1m2, 0.0);
	Vector<double> O1(m1m2 * m1m2, 0.0);

	Vector<int> QN_up(m1m2, 0), QN_down(m1m2, 0);
	Vector<long> Multi_Index(m1m2);

	CSparseSymmetric<double, long> SubHSB;
	CSparseSymmetric<double, long> H12;
	CSparseUnitary<double, long> __c1a;
	CSparseUnitary<double, long> __c1b;

	FILE* stream = fopen("DAT\\InfDMRG_Hub_QN.dat", "w+");
		fprintf(stream, "DMRG RESULTS:");
	fclose(stream);

	int Step = -1;

	while (fabs(Convergence) > m_EPS_DMRG)
	{
		++Step;
		L = SitesNumber;
		m1m2 = m1 * m2;
		m_Size = m1 * m2 * m2 * m1;

		int N_up;
		int N_down;

		if ((qn_up == 2 && qn_down == 2) || (qn_up == 4 && qn_down == 0))
		{
			N_up = SitesNumber * qn_up / 4;
			N_down = SitesNumber * qn_down / 4;
		}
		else
		{
			N_up = SitesNumber - (4 - qn_up);
			N_down = qn_down;
		}

		int SubSize = GetSubSize(mi, N_up, N_down);

		if (SubSize > MaxSize)
		{
			MaxSize = SubSize;
			TargetState.Realloc(SubSize);
		}
			
		SubHSB.reserve(SubSize, SubSize, SubSize);
		SubHSB.resize(SubSize, 0, 0);

		cout << "\nIteration = " << (Step + 1) << endl;
		cout << "N_up = " << N_up << "; N_down = " << N_down << endl;
		cout << "mi : " << mi[0] << ", " << mi[1] << ", " << mi[2] << ", " << mi[3] << endl;
		cout << SitesNumber << " sites" << endl;
		cout << "Size = " << m_Size;
		cout << ";\t SubSize = " << SubSize << endl;

		FormSuperBlock(mi, SubHSB, Step);
		
		cout << "Sub : NEL = " << SubHSB.ir.size();
		cout << ";\t Sparsity = " << SubHSB.GetSparsity() << " %" << endl;

		GSE = Diagonalize_HSB(SubSize, SubHSB, TargetState);

		Convergence = GSE_prev - GSE / SitesNumber;
		GSE_prev = GSE / SitesNumber;

		cout << "Sub : GSE = " << GSE;
		cout << "; GSE / L = " << (GSE / SitesNumber);
		cout << "; Convergence = " << Convergence << endl;

		stream = fopen("DAT\\InfDMRG_Hub_QN.dat", "a+");
			fprintf(stream, "\nE(L=%d,\tm=%d)\t= %lf;\tAccuracy = %6.5e;\tIteration = %d",
					SitesNumber, mi[0], (GSE / SitesNumber), m_Accuracy, Step + 1);
		fclose(stream);

//		if (_getche() == 3)	exit(0);
		int _m1 = (m1m2 > m_m) ? m_m : m1m2;

		FormDensityMatrix12(mi, DM12, TargetState, QN_up, QN_down);
		Arrange_ro(m1m2, QN_up, QN_down, *QN_up_1, *QN_down_1, DM12, Multi_Index);

		m_Accuracy = Diagonalize_ro(m1m2, DM12, w1, O1, QN_up, QN_down, Multi_Index);
		Form12BlockSystem(mi, H12, __c1a, __c1b, Multi_Index, Step);
		Arrange_QN(_m1, m1m2, O1, QN_up, QN_down, *QN_up_1, *QN_down_1);

		Truncation(_m1, m1m2, DM12, H12, O1, *H1);
		Truncation(_m1, m1m2, DM12, __c1a, O1, *c1a);
		Truncation(_m1, m1m2, DM12, __c1b, O1, *c1b);


		FormDensityMatrix34(mi, DM12, TargetState, QN_up, QN_down);
		Arrange_ro(m1m2, QN_up, QN_down, *QN_up_4, *QN_down_4, DM12, Multi_Index);

		m_Accuracy = Diagonalize_ro(m1m2, DM12, w1, O1, QN_up, QN_down, Multi_Index);
		Form34BlockSystem(mi, H12, __c1a, __c1b, Multi_Index, Step);
		Arrange_QN(_m1, m1m2, O1, QN_up, QN_down, *QN_up_4, *QN_down_4);

		Truncation(_m1, m1m2, DM12, H12, O1, *H4);
		Truncation(_m1, m1m2, DM12, __c1a, O1, *c4a);
		Truncation(_m1, m1m2, DM12, __c1b, O1, *c4b);

		m1 = (m1m2 > m_m) ? m_m : m1m2;
		mi[0] = mi[3] = m1;

		SitesNumber += 2;
	}

	return GSE;
}

void CInfDMRG_Hub_QN::Initialization(double U)
{
	for (int i = 0; i < N; i++)
	{
		(*QN_up_1)[i] = nna[i];
		  QN_up_2 [i] = nna[i];
		  QN_up_3 [i] = nna[i];
		(*QN_up_4)[i] = nna[i];

		(*QN_down_1)[i] = nnb[i];
		  QN_down_2 [i] = nnb[i];
		  QN_down_3 [i] = nnb[i];
		(*QN_down_4)[i] = nnb[i];

		(*H1)[i][i] = U * (nna[i] * nnb[i]);
		  H2 [i][i] = U * (nna[i] * nnb[i]);
		  H3 [i][i] = U * (nna[i] * nnb[i]);
		(*H4)[i][i] = U * (nna[i] * nnb[i]);

		for (int j = 0; j < N; j++)
		{
			(*c1a)[i][j] = __cna[i][j];
			  c2a [i][j] = cna[i][j];
			  c3a [i][j] = __cna[i][j];
			(*c4a)[i][j] = cna[i][j];

			(*c1b)[i][j] = __cnb[i][j];
			  c2b [i][j] = cnb[i][j];
			  c3b [i][j] = __cnb[i][j];
			(*c4b)[i][j] = cnb[i][j];
		}
	}
}

int CInfDMRG_Hub_QN::GetSubSize(int* mi, int qn_up, int qn_down)
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
		int L = k1 * _m234 + k2 * _m34 + k3 * _m4 + k4;
		
		int QN_up_H = (*QN_up_1)[k1] + QN_up_2[k2] + QN_up_3[k3] + (*QN_up_4)[k4];
		int QN_down_H = (*QN_down_1)[k1] + QN_down_2[k2] + QN_down_3[k3] + (*QN_down_4)[k4];

		if (QN_up_H == qn_up && QN_down_H == qn_down)
		{
			QN_index_H.push_back(L);
			Counter++;
		}
	}

	return Counter;
}

void CInfDMRG_Hub_QN::FormSuperBlock(int* mi, CSparseSymmetric<double, long>& SubHSB, int Step)
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
				Summa_t += (*c1a)[j1][i1] * c2a[i2][j2] + c2a[j2][i2] * (*c1a)[i1][j1];
				Summa_t += (*c1b)[j1][i1] * c2b[i2][j2] + c2b[j2][i2] * (*c1b)[i1][j1];
			}

			if (i1 == j1 && i4 == j4)
			{
				Summa_t += __cna[j2][i2] * cna[i3][j3] + cna[j3][i3] * __cna[i2][j2];
				Summa_t += __cnb[j2][i2] * cnb[i3][j3] + cnb[j3][i3] * __cnb[i2][j2];
			}

			if (i1 == j1 && i2 == j2)
			{
				Summa_t += c3a[j3][i3] * (*c4a)[i4][j4] + (*c4a)[j4][i4] * c3a[i3][j3];
				Summa_t += c3b[j3][i3] * (*c4b)[i4][j4] + (*c4b)[j4][i4] * c3b[i3][j3];
			}

			double Summa = Summa_t * m_t1 + Summa_B;
			
			SubHSB.push_back(Summa, Index_L, Index_R);
		}
	}
}

double CInfDMRG_Hub_QN::Diagonalize_HSB(int Size, CSparseSymmetric<double, long>& HSB, Vector<double>& TargetState)
{
	double E = 0.0;

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
		int		nDim = 4 * NEVEC;

		Vector<double>	EVAL(NEVEC);

		Davidson(NEVEC, SA, EVAL, TargetState, nDim, m_EPS_Diag);

		E = EVAL[0];
	}

	return E;
}

void CInfDMRG_Hub_QN::FormDensityMatrix12(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Vector<int>& QN_up, Vector<int>& QN_down)
{
	int _m234 = mi[1] * mi[2] * mi[3];
	int _m34 = mi[2] * mi[3];
	int _m4 = mi[3];

	int SubSize = QN_index_H.size();
	int m1m2 = mi[0] * mi[1];

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

	for (int i = 0; i < m1m2; i++)
	{
		DM[i * m1m2 + i] = 0.0;

		int i1 = i / mi[1];
		int i2 = i % mi[1];

		QN_up  [i] = (*QN_up_1)  [i1] + QN_up_2  [i2];
		QN_down[i] = (*QN_down_1)[i1] + QN_down_2[i2];

		for (int j = i + 1; j < m1m2; j++)
		{
			DM[j * m1m2 + i] = DM[i * m1m2 + j] = 0.0;
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

			DM[RR * m1m2 + LL] = DM[LL * m1m2 + RR] = Summa;

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

void CInfDMRG_Hub_QN::FormDensityMatrix34(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Vector<int>& QN_up, Vector<int>& QN_down)
{
	int _m234 = mi[1] * mi[2] * mi[3];
	int _m34 = mi[2] * mi[3];
	int _m4 = mi[3];

	int SubSize = QN_index_H.size();
	int m3m4 = mi[2] * mi[3];

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

	for (int i = 0; i < m3m4; i++)
	{
		DM[i * m3m4 + i] = 0.0;

		int i3 = i / mi[3];
		int i4 = i % mi[3];

		QN_up  [i] = QN_up_3  [i3] + (*QN_up_4)  [i4];
		QN_down[i] = QN_down_3[i3] + (*QN_down_4)[i4];

		for (int j = i + 1; j < m3m4; j++)
		{
			DM[j * m3m4 + i] = DM[i * m3m4 + j] = 0.0;
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

			DM[RR * m3m4 + LL] = DM[LL * m3m4 + RR] = Summa;

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

void Orthogonality(int N, int LDA, int NEVEC, double EPS, Vector<double>& EVAL,
				   Vector<double>& EVEC)
{
	for (int j = 0; j < NEVEC - 1; j++)
	{
		for (int k = j + 1; k < NEVEC; k++)
		{
			double Summa = 0.0;

			for (int i = 0; i < N; i++)
			{
				Summa += EVEC[j * LDA + i] * EVEC[k * LDA + i];
			}

			if (fabs(Summa) > EPS)
			{
				cout << "Non-orthogonal vectors : ";
				cout << j << " * " << k << " = " << Summa << endl;
//				if (_getche() == 3)	exit(0);
			}
		}
	}
}

double CInfDMRG_Hub_QN::Diagonalize_ro(int Size, Vector<double>& DM, Vector<double>& w, Vector<double>& O,
					Vector<int>& QN_up, Vector<int>& QN_down, Vector<long>& Multi_Index)
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

		if (QN_up[i] != QN_up[i + 1] || QN_down[i] != QN_down[i + 1])
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

	sort_qn_2_c(Size, SMALL, w.getVec(), O.getVec(), QN_up.getVec(), QN_down.getVec(), Multi_Index.getVec(), index.getVec());

	double Accuracy = 0.0;
	double Spur = 0.0;
	double Spur_O = 0.0;

	for (int i = 0; i < __m; i++)
	{
		Accuracy += w[i];
	}

	for (int i = 0; i < Size; i++)
	{
		Spur += w[i];

		for (int k = 0; k < Size; k++)
		{
			Spur_O += O[i * Size + k] * O[i * Size + k];
		}
	}

	Accuracy = Spur - Accuracy;

	cout << "Spur(O) = " << Spur_O;
	cout << ";\t Spur(ro) = " << Spur << ";\t Accuracy = " << Accuracy << endl;

	if (fabs(Spur_O - Size) > 1.0E-008 || fabs(Spur - 1) > 1.0E-008)
	{
		cout << "Spur(O) - Size = " << (Spur_O - Size);
		cout << ";\t Spur(ro) - 1 = " << (Spur - 1) << endl;
		if (_getche() == 3)	exit(0);
	}

	return Accuracy;
}

void CInfDMRG_Hub_QN::Form12BlockSystem(int* mi, CSparseSymmetric<double, long>& H12,
					CSparseUnitary<double, long>& __c1a,
					CSparseUnitary<double, long>& __c1b,
					Vector<long>& Multi_Index, int Step)
{
	int m1m2 = mi[0] * mi[1];

	H12.reserve(m1m2, m1m2, m1m2);
	H12.resize(m1m2, 0, 0);

	__c1a.resize(m1m2 + 1, m1m2 + 1);
	__c1a.ija[0] = m1m2 + 1;

	__c1b.resize(m1m2 + 1, m1m2 + 1);
	__c1b.ija[0] = m1m2 + 1;

	for (int i = 0; i < m1m2; i++)
	{
		long L = Multi_Index[i];
		int i1 = L / mi[1];
		int i2 = L % mi[1];
	
		__c1a.sa[i] = 0.0;
		__c1b.sa[i] = 0.0;

		for (int j = 0; j < m1m2; j++)
		{
			long R = Multi_Index[j];
			int j1 = R / mi[1];
			int j2 = R % mi[1];
		
			if (i1 == j1)
			{
				__c1a.push_back(__cna[i2][j2], i, j);
				__c1b.push_back(__cnb[i2][j2], i, j);
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

			Summa_t += (*c1a)[j1][i1] * c2a[i2][j2] + c2a[j2][i2] * (*c1a)[i1][j1];
			Summa_t += (*c1b)[j1][i1] * c2b[i2][j2] + c2b[j2][i2] * (*c1b)[i1][j1];

			double Summa = Summa_t * m_t1 + Summa_B;

			H12.push_back(Summa, i, j);
		}

		__c1a.ija[i + 1] = __c1a.ija.size();
		__c1b.ija[i + 1] = __c1b.ija.size();
	}
}

void CInfDMRG_Hub_QN::Form34BlockSystem(int* mi, CSparseSymmetric<double, long>& H34,
					CSparseUnitary<double, long>& __c4a,
					CSparseUnitary<double, long>& __c4b,
					Vector<long>& Multi_Index, int Step)
{
	int m3m4 = mi[2] * mi[3];

	H34.reserve(m3m4, m3m4, m3m4);
	H34.resize(m3m4, 0, 0);

	__c4a.resize(m3m4 + 1, m3m4 + 1);
	__c4a.ija[0] = m3m4 + 1;

	__c4b.resize(m3m4 + 1, m3m4 + 1);
	__c4b.ija[0] = m3m4 + 1;

	for (int i = 0; i < m3m4; i++)
	{
		long L = Multi_Index[i];
		int i3 = L / mi[3];
		int i4 = L % mi[3];

		__c4a.sa[i] = 0.0;
		__c4b.sa[i] = 0.0;

		for (int j = 0; j < m3m4; j++)
		{
			long R = Multi_Index[j];
			int j3 = R / mi[3];
			int j4 = R % mi[3];
		
			if (i4 == j4)
			{
				__c4a.push_back(cna[i3][j3], i, j);
				__c4b.push_back(cnb[i3][j3], i, j);
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

			Summa_t += c3a[j3][i3] * (*c4a)[i4][j4] + (*c4a)[j4][i4] * c3a[i3][j3];
			Summa_t += c3b[j3][i3] * (*c4b)[i4][j4] + (*c4b)[j4][i4] * c3b[i3][j3];

			double Summa = Summa_t * m_t1 + Summa_B;

			H34.push_back(Summa, i, j);
		}

		__c4a.ija[i + 1] = __c4a.ija.size();
		__c4b.ija[i + 1] = __c4b.ija.size();
	}
}

void CInfDMRG_Hub_QN::Arrange_ro(int m1m2, Vector<int>& Source_up, Vector<int>& Source_down,
				Vector<int>& Destin_up, Vector<int>& Destin_down,
				Vector<double>& DM, Vector<long>& Multi_Index)
{
	for (int i = 0; i < m1m2; i++)
	{
		Destin_up[i] = Source_up[i];
		Destin_down[i] = Source_down[i];
		Multi_Index[i] = i;
	}

//	cout << "\nDM before horizontal arrangement:" << endl;
//	Debug(m1m2, m1m2, DM, 0.0);
//	cout << endl;
//	if (_getche() == 3)	exit(0);

	for (int i = 0; i < m1m2 - 1; i++)
	{
		int k = i;

		double qn_up = Source_up[i];
		double qn_down = Source_down[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (Source_up[j] > qn_up ||
				(Source_up[j] == qn_up && Source_down[j] > qn_down))
			{
				k = j;
				qn_up = Source_up[j];
				qn_down = Source_down[j];
			}
		}

		if (k != i)
		{
			Source_up[k] = Source_up[i];
			Source_up[i] = qn_up;

			Source_down[k] = Source_down[i];
			Source_down[i] = qn_down;

			long L = Multi_Index[i];
			Multi_Index[i] = Multi_Index[k];
			Multi_Index[k] = L;

			for (int j = 0; j < m1m2; j++)
			{
				qn_up = DM[j * m1m2 + i];
				DM[j * m1m2 + i] = DM[j * m1m2 + k];
				DM[j * m1m2 + k] = qn_up;
			}
		}
	}

//	cout << "Multi_Index :" << endl;
//	Debug(m1m2, Multi_Index);
//	if (_getche() == 3)	exit(0);

//	cout << "\nDM after horizontal and before vertical arrangement:" << endl;
//	Debug(m1m2, m1m2, DM, 0.0);
//	cout << endl;
//	if (_getche() == 3)	exit(0);

	for (int i = 0; i < m1m2 - 1; i++)
	{
		int k = i;

		double qn_up = Destin_up[i];
		double qn_down = Destin_down[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (Destin_up[j] > qn_up ||
				(Destin_up[j] == qn_up && Destin_down[j] > qn_down))
			{
				k = j;
				qn_up = Destin_up[j];
				qn_down = Destin_down[j];
			}
		}

		if (k != i)
		{
			Destin_up[k] = Destin_up[i];
			Destin_up[i] = qn_up;

			Destin_down[k] = Destin_down[i];
			Destin_down[i] = qn_down;

			for (int j = 0; j < m1m2; j++)
			{
				qn_up = DM[i * m1m2 + j];
				DM[i * m1m2 + j] = DM[k * m1m2 + j];
				DM[k * m1m2 + j] = qn_up;
			}
		}
	}

//	cout << "\nDM after vertical arrangement:" << endl;
//	Debug(m1m2, m1m2, DM, 0.0);
//	cout << endl;
//	if (_getche() == 3)	exit(0);
}

void CInfDMRG_Hub_QN::Arrange_QN(int m, int m1m2, Vector<double>& O,
				Vector<int>& Source_up, Vector<int>& Source_down,
				Vector<int>& Destin_up, Vector<int>& Destin_down)
{
	for (int i = 0; i < m; i++)
	{
		Destin_up  [i] = Source_up  [i];
		Destin_down[i] = Source_down[i];
	}
}

void CInfDMRG_Hub_QN::Truncation(int m, int m1m2, Vector<double>& H12_O, CSparse<double, long>& H12,
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
