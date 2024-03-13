//
// File InfDMRG_S12_QN.cpp
//
// Implementation of class CInfDMRG_S12_QN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Heisenberg S=1/2 model
//
#include <iostream>
#include <vector>
#include <conio.h>
#include <time.h>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"

#include "Spin2.h"
//#include "Spin3.h"

#include "SpinOp.h"
#include "CSparse.hpp"
#include "InfDMRG_S12_QN.h"

#include "Debug.hpp"

#if __SPIN == 2
	const int N = 2;
#elif __SPIN == 3
	const int N = 3;
#endif

Matrix<double> S_up_2(N, N), S_up_3(N, N);
Matrix<double> S_z_2(N, N), S_z_3(N, N);
Matrix<double> H2(N, N), H3(N, N);
Vector<double> QN2(N), QN3(N);

vector<long, allocator<long> > QN_index_H;

vector<long, allocator<long> > vec_1;
vector<long, allocator<long> > vec_2;
vector<long, allocator<long> > vec_3;
vector<long, allocator<long> > vec_4;


CInfDMRG_S12_QN::CInfDMRG_S12_QN(int m, double J1, double J2, double S)
{
	if (N != (int)(2 * S + 1)) Error("Wrong spin");

	m_m = m;
	m_J1 = J1;
	m_J2 = J2;
	m_Size = m_m * N * N * m_m;
	m_Accuracy = 0.0;
	m_EPS_DMRG = 1.0E-005;
	m_EPS_Diag = 1.0E-008;
	m_Iterations = 100;

	S_up_1 = new Matrix<double>(m_m, m_m, 0.0);
	S_up_4 = new Matrix<double>(m_m, m_m, 0.0);
	S_z_1 = new Matrix<double>(m_m, m_m, 0.0);
	S_z_4 = new Matrix<double>(m_m, m_m, 0.0);
	
	H1 = new Matrix<double>(m_m, m_m, 0.0);
	H4 = new Matrix<double>(m_m, m_m, 0.0);
	
	QN1 = new Vector<double>(m_m * N, 0.0);
	QN4 = new Vector<double>(m_m * N, 0.0);
	
	Initialization();
};

CInfDMRG_S12_QN::~CInfDMRG_S12_QN()
{
	delete QN4;
	delete QN1;
	
	delete H4;
	delete H1;
	
	delete S_z_4;
	delete S_z_1;
	delete S_up_4;
	delete S_up_1;
}

int CInfDMRG_S12_QN::IsExit(int Param, double Convergence, int Iteration)
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

double CInfDMRG_S12_QN::Iterations(int& L, Vector<double>& TargetState, double qn_Sz, int Param)
{
	int m1 = N, m2 = N;
	int m1m2 = m_m * m2;
	int mi[4] = { m1, m2, m2, m1 };
	int SitesNumber = 4;

	double GSE, GSE_per_site = 0.0, Convergence = 1.0;
	long MaxSize = TargetState.getSize();

	Vector<double> DM12(m1m2 * m1m2, 0.0);
	Vector<double> w1(m1m2, 0.0);
	Vector<double> O1(m1m2 * m1m2, 0.0);

	Vector<double>	QN(m1m2, 0.0);
	Vector<long>	Multi_Index(m1m2);

	CSparseSymmetric<double, long> SubHSB;
	CSparseSymmetric<double, long> H12;
	CSparseUnitary<double, long> _S_up_1;
	CSparseUnitary<double, long> _S_z_1;
	
	struct tm *newtime;
	time_t start, finish;

	time(&start);
	newtime = localtime(&start);

	FILE* stream = fopen("DAT\\InfDMRG_S12_QN.dat", "w+");
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

#if __SPIN == 2
		double Sz_total = (qn_Sz == 0 || qn_Sz == 2) ? SitesNumber * qn_Sz / 4 : qn_Sz;
#elif __SPIN == 3
		double Sz_total = (qn_Sz == 0 || qn_Sz == 4) ? SitesNumber * qn_Sz / 4 : qn_Sz;
#endif

		int SubSize = GetSubSize(mi, Sz_total);
		
		if (SubSize > MaxSize)
		{
			MaxSize = SubSize;
			TargetState.Realloc(SubSize);
		}
			
		SubHSB.reserve(SubSize, SubSize, SubSize);
		SubHSB.resize(SubSize, 0, 0);

		cout << "\nIteration = " << Iteration << endl;
		cout << "Sz_total = " << Sz_total << endl;
		cout << "mi : " << mi[0] << ", " << mi[1] << ", " << mi[2] << ", " << mi[3] << endl;
		cout << SitesNumber << " sites" << endl;
		cout << "Size = " << m_Size;
		cout << ";\t SubSize = " << SubSize << endl;
//		if (_getche() == 3)	exit(0);

		FormSuperBlock(mi, SubHSB, Iteration);

		cout << "Sub : NEL = " << SubHSB.ir.size();
		cout << ";\t Sparsity = " << SubHSB.GetSparsity() << " %" << endl;

		GSE = Diagonalize_HSB(SubSize, SubHSB, TargetState);

		Convergence = GSE_per_site - GSE / SitesNumber;
		GSE_per_site = GSE / SitesNumber;

		cout << "Sub : GSE = " << GSE;
		cout << "; GSE / L = " << GSE_per_site;
		cout << "; Convergence = " << Convergence << endl;

		stream = fopen("DAT\\InfDMRG_S12_QN.dat", "a+");
			fprintf(stream, "\nE(L=%d,\tm=%d)", SitesNumber, mi[0]);
			fprintf(stream, "\t= %lf;", GSE);
			fprintf(stream, "\tE / L = %lf;", GSE_per_site);
			fprintf(stream, "\tTruncation error = %6.5e;", m_Accuracy);
			fprintf(stream, "\tIteration = %d;", Iteration);
			fprintf(stream, "\tConvergence = %6.5e;", Convergence);
		fclose(stream);

		int _m1 = (m1m2 > m_m) ? m_m : m1m2;

//////////////////////////////////////////////////////////////////////////////////
		
		FormDensityMatrix12(mi, DM12, TargetState, QN);
		Arrange_ro(m1m2, QN, *QN1, DM12, Multi_Index);
		
		m_Accuracy = Diagonalize_ro(m1m2, DM12, w1, O1, QN, Multi_Index);
		Form12BlockSystem(mi, H12, _S_up_1, _S_z_1, Multi_Index, Iteration);
		Arrange_QN(_m1, m1m2, O1, QN, *QN1);

		Vector<double> One(m1m2 * m1m2, 0.0);
		for (int i = 0; i < m1m2; i++) One[i * m1m2 + i] = 1.0;
		H12.Multiply_TR(m1m2, _m1, DM12, One);

		Truncation(_m1, m1m2, DM12, H12, O1, *H1);
		Truncation(_m1, m1m2, DM12, _S_up_1, O1, *S_up_1);
		Truncation(_m1, m1m2, DM12, _S_z_1, O1, *S_z_1);

//////////////////////////////////////////////////////////////////////////////////

		FormDensityMatrix34(mi, DM12, TargetState, QN);
		Arrange_ro(m1m2, QN, *QN4, DM12, Multi_Index);
		
		m_Accuracy = Diagonalize_ro(m1m2, DM12, w1, O1, QN, Multi_Index);
		Form34BlockSystem(mi, H12, _S_up_1, _S_z_1, Multi_Index, Iteration);
		Arrange_QN(_m1, m1m2, O1, QN, *QN4);

		Truncation(_m1, m1m2, DM12, H12, O1, *H4);
		Truncation(_m1, m1m2, DM12, _S_up_1, O1, *S_up_4);
		Truncation(_m1, m1m2, DM12, _S_z_1, O1, *S_z_4);

//////////////////////////////////////////////////////////////////////////////////

		m1 = (m1m2 > m_m) ? m_m : m1m2;
		mi[0] = mi[3] = m1;

		SitesNumber += 2;
	}

	time(&finish);
	newtime = localtime(&finish);
	double elapsed_time = difftime(finish, start);
	
	int hours = (int)(elapsed_time / 3600.0);
	int minutes = (int)((elapsed_time - hours * 3600) / 60.0);
	int seconds = (int)(elapsed_time - hours * 3600 - minutes * 60);

	stream = fopen("DAT\\InfDMRG_S12_QN.dat", "a+");
		fprintf(stream, "\n\nFINISH: %s", asctime(newtime));
		fprintf(stream, "\nELAPSED TIME: %d hours %d minutes %d seconds",
				hours, minutes, seconds);
	fclose(stream);

	return GSE;
}

void CInfDMRG_S12_QN::Initialization(void)
{
	for (int i = 0; i < N; i++)
	{
		(*QN1)[i] = Sn_z[i][i];
		  QN2 [i] = Sn_z[i][i];
		  QN3 [i] = Sn_z[i][i];
		(*QN4)[i] = Sn_z[i][i];

		for (int j = 0; j < N; j++)
		{
			(*S_up_1)[i][j] = Sn_up[i][j];
			  S_up_2 [i][j] = Sn_up[i][j];
			  S_up_3 [i][j] = Sn_up[i][j];
			(*S_up_4)[i][j] = Sn_up[i][j];

			(*S_z_1)[i][j] = Sn_z[i][j];
			  S_z_2 [i][j] = Sn_z[i][j];
			  S_z_3 [i][j] = Sn_z[i][j];
			(*S_z_4)[i][j] = Sn_z[i][j];

			(*H1)[i][j] = 0.0;
			  H2 [i][j] = 0.0;
			  H3 [i][j] = 0.0;
			(*H4)[i][j] = 0.0;
		}
	}
}

int CInfDMRG_S12_QN::GetSubSize(int* mi, double qn)
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
		
		double QN = (*QN1)[k1] + QN2[k2] + QN3[k3] + (*QN4)[k4];
		
		if (QN == qn)
		{
			QN_index_H.push_back(L);
			Counter++;
		}
	}

	return Counter;
}

void CInfDMRG_S12_QN::FormSuperBlock(int* mi, CSparseSymmetric<double, long>& SubHSB, int Iteration)
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

			double Summa_J1 = 0.0;
			double Summa_J2 = 0.0;

			if (i3 == j3 && i4 == j4)
			{
				Summa_J1 += (*S_z_1)[i1][j1] * S_z_2[i2][j2];
				Summa_J1 += 0.5 * (*S_up_1)[i1][j1] * S_up_2[j2][i2];
				Summa_J1 += 0.5 * (*S_up_1)[j1][i1] * S_up_2[i2][j2];
			}

			if (i1 == j1 && i4 == j4)
			{
				Summa_J2 += S_z_2[i2][j2] * S_z_3[i3][j3];
				Summa_J2 += 0.5 * S_up_2[i2][j2] * S_up_3[j3][i3];
				Summa_J2 += 0.5 * S_up_2[j2][i2] * S_up_3[i3][j3];
			}

			if (i1 == j1 && i2 == j2)
			{
				Summa_J1 += S_z_3[i3][j3] * (*S_z_4)[i4][j4];
				Summa_J1 += 0.5 * S_up_3[i3][j3] * (*S_up_4)[j4][i4];
				Summa_J1 += 0.5 * S_up_3[j3][i3] * (*S_up_4)[i4][j4];
			}

			double Summa = 0.0;

			if (Iteration % 2 == 0)
			{
				Summa = Summa_J1 * m_J2 + Summa_J2 * m_J1 + Summa_B;
			}
			else
			{
				Summa = Summa_J1 * m_J1 + Summa_J2 * m_J2 + Summa_B;
			}

			SubHSB.push_back(Summa, Index_L, Index_R);
		}
	}
}

double CInfDMRG_S12_QN::Diagonalize_HSB(int Size, CSparseSymmetric<double, long>& HSB, Vector<double>& TargetState)
{
	double E;

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

		int	NEVEC = 1;
		int	nDim = 25 * NEVEC;

		Vector<double>	EVAL(NEVEC);

		Davidson(NEVEC, SA, EVAL, TargetState, nDim, m_EPS_Diag);

		E = EVAL[0];
	}

	return E;
}

void CInfDMRG_S12_QN::FormDensityMatrix12(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Vector<double>& QN)
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

		QN[i] = (*QN1)[i1] + QN2[i2];

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

void CInfDMRG_S12_QN::FormDensityMatrix34(int* mi, Vector<double>& DM, Vector<double>& TargetState,
					Vector<double>& QN)
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

		QN[i] = QN3[i3] + (*QN4)[i4];

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

double CInfDMRG_S12_QN::Diagonalize_ro(int Size, Vector<double>& DM, Vector<double>& w, Vector<double>& O,
					Vector<double>& QN, Vector<long>& Multi_Index)
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

		if (QN[i] != QN[i + 1])
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

	sort_qn_c(Size, SMALL, w.getVec(), O.getVec(), QN.getVec(), Multi_Index.getVec(), index.getVec());

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

void CInfDMRG_S12_QN::Form12BlockSystem(int* mi, CSparseSymmetric<double, long>& H12,
					CSparseUnitary<double, long>& _S_up_1,
					CSparseUnitary<double, long>& _S_z_1,
					Vector<long>& Multi_Index, int Iteration)
{
	int m1m2 = mi[0] * mi[1];

	H12.reserve(m1m2, m1m2, m1m2);
	H12.resize(m1m2, 0, 0);

	_S_up_1.resize(m1m2 + 1, m1m2 + 1);
	_S_up_1.ija[0] = m1m2 + 1;

	_S_z_1.resize(m1m2 + 1, m1m2 + 1);
	_S_z_1.ija[0] = m1m2 + 1;

	for (int i = 0; i < m1m2; i++)
	{
		long L = Multi_Index[i];
		int i1 = L / mi[1];
		int i2 = L % mi[1];
	
		_S_up_1.sa[i] = 0.0;
		_S_z_1.sa[i] = 0.0;

		for (int j = 0; j < m1m2; j++)
		{
			long R = Multi_Index[j];
			int j1 = R / mi[1];
			int j2 = R % mi[1];
		
			if (i1 == j1)
			{
				_S_up_1.push_back(S_up_2[i2][j2], i, j);
				_S_z_1.push_back(S_z_2[i2][j2], i, j);
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

			double Summa_J = 0.0;

			Summa_J += (*S_z_1)[i1][j1] * S_z_2[i2][j2];
			Summa_J += 0.5 * (*S_up_1)[i1][j1] * S_up_2[j2][i2];
			Summa_J += 0.5 * (*S_up_1)[j1][i1] * S_up_2[i2][j2];

			double Summa = 0.0;

			if (Iteration % 2 == 0)
			{
				Summa = Summa_J * m_J2 + Summa_B;
			}
			else
			{
				Summa = Summa_J * m_J1 + Summa_B;
			}

			H12.push_back(Summa, i, j);
		}

		_S_up_1.ija[i + 1] = _S_up_1.ija.size();
		_S_z_1.ija[i + 1] = _S_z_1.ija.size();
	}
}

void CInfDMRG_S12_QN::Form34BlockSystem(int* mi, CSparseSymmetric<double, long>& H34,
					CSparseUnitary<double, long>& _S_up_4,
					CSparseUnitary<double, long>& _S_z_4,
					Vector<long>& Multi_Index, int Iteration)
{
	int m3m4 = mi[2] * mi[3];

	H34.reserve(m3m4, m3m4, m3m4);
	H34.resize(m3m4, 0, 0);

	_S_up_4.resize(m3m4 + 1, m3m4 + 1);
	_S_up_4.ija[0] = m3m4 + 1;

	_S_z_4.resize(m3m4 + 1, m3m4 + 1);
	_S_z_4.ija[0] = m3m4 + 1;

	for (int i = 0; i < m3m4; i++)
	{
		long L = Multi_Index[i];
		int i3 = L / mi[3];
		int i4 = L % mi[3];

		_S_up_4.sa[i] = 0.0;
		_S_z_4.sa[i] = 0.0;

		for (int j = 0; j < m3m4; j++)
		{
			long R = Multi_Index[j];
			int j3 = R / mi[3];
			int j4 = R % mi[3];
		
			if (i4 == j4)
			{
				_S_up_4.push_back(S_up_3[i3][j3], i, j);
				_S_z_4.push_back(S_z_3[i3][j3], i, j);
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

			double Summa_J = 0.0;

			Summa_J += S_z_3[i3][j3] * (*S_z_4)[i4][j4];
			Summa_J += 0.5 * S_up_3[i3][j3] * (*S_up_4)[j4][i4];
			Summa_J += 0.5 * S_up_3[j3][i3] * (*S_up_4)[i4][j4];

			double Summa = 0.0;

			if (Iteration % 2 == 0)
			{
				Summa = Summa_J * m_J2 + Summa_B;
			}
			else
			{
				Summa = Summa_J * m_J1 + Summa_B;
			}

			H34.push_back(Summa, i, j);
		}

		_S_up_4.ija[i + 1] = _S_up_4.ija.size();
		_S_z_4.ija[i + 1] = _S_z_4.ija.size();
	}
}

void CInfDMRG_S12_QN::Arrange_ro(int m1m2, Vector<double>& QN_Source, Vector<double>& QN_Destin,
				Vector<double>& DM, Vector<long>& Multi_Index)
{
	for (int i = 0; i < m1m2; i++)
	{
		QN_Destin[i] = QN_Source[i];
		Multi_Index[i] = i;
	}

//	cout << "\nDM before horizontal arrangement:" << endl;
//	Debug(m1m2, m1m2, DM, 0.0);
//	cout << endl;
//	if (_getche() == 3)	exit(0);

	for (int i = 0; i < m1m2 - 1; i++)
	{
		int k = i;

		double qn = QN_Source[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (QN_Source[j] > qn)
			{
				k = j;
				qn = QN_Source[j];
			}
		}

		if (k != i)
		{
			QN_Source[k] = QN_Source[i];
			QN_Source[i] = qn;

			long L = Multi_Index[i];
			Multi_Index[i] = Multi_Index[k];
			Multi_Index[k] = L;

			for (int j = 0; j < m1m2; j++)
			{
				qn = DM[j * m1m2 + i];
				DM[j * m1m2 + i] = DM[j * m1m2 + k];
				DM[j * m1m2 + k] = qn;
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

		double qn = QN_Destin[i];

		for (int j = i + 1; j < m1m2; j++)
		{
			if (QN_Destin[j] > qn)
			{
				k = j;
				qn = QN_Destin[j];
			}
		}

		if (k != i)
		{
			QN_Destin[k] = QN_Destin[i];
			QN_Destin[i] = qn;

			for (int j = 0; j < m1m2; j++)
			{
				qn = DM[i * m1m2 + j];
				DM[i * m1m2 + j] = DM[k * m1m2 + j];
				DM[k * m1m2 + j] = qn;
			}
		}
	}

//	cout << "\nDM after vertical arrangement:" << endl;
//	Debug(m1m2, m1m2, DM, 0.0);
//	cout << endl;
//	if (_getche() == 3)	exit(0);
}

void CInfDMRG_S12_QN::Arrange_QN(int m, int m1m2, Vector<double>& O, Vector<double>& Source, Vector<double>& Destin)
{
	for (int i = 0; i < m; i++)
	{
		Destin[i] = Source[i];
	}
}

void CInfDMRG_S12_QN::Truncation(int m, int m1m2, Vector<double>& H12_O, CSparse<double, long>& H12,
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
