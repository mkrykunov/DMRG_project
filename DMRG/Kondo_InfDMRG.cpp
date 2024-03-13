// g++ -o Kondo_InfDMRG.exe Kondo_InfDMRG.cpp External.cpp InfDMRGKondoQN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_3_c.c Davidson.cpp
//
#include <iostream>
#include <vector>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "CSparse.hpp"
#include "InfDMRGKondoQN.h"


int main(void)
{
	int L = 20;
	int m = 32;
	int N = 32;
	int Param = 0;
	int Krylov = 25;

	double t = -1.0;
	double alpha = 10 * t;
	double UL = 10 * t;
	double UA = 10 * t;
	double J = -39.0;
	
	int qn_LA = 1.0;
	int qn_d = 1.0;
	double qn_S = 4.0;

	FILE* stream;

	if ((stream = fopen("DAT\\Kondo.dat", "r")) != NULL)
	{
		fscanf(stream, "%lf %lf %lf %lf %lf %d %d %lf %d %d",
			   &t, &alpha, &UL, &UA, &J, &qn_LA, &qn_d, &qn_S, &Param, &Krylov);
		fclose(stream);
	}

	t /= alpha;
	UL /= alpha;
	UA /= alpha;
	J /= fabs(alpha);
	alpha /= alpha;

	cout << "t = " << t << ";   alpha = " << alpha << ";   UL = " << UL;
	cout << ";   UA = " << UA << ";   J = " << J << endl;
	cout << "qn_LA = " << qn_LA << ";   qn_d = " << qn_d;
	cout << ";   qn_S = " << qn_S << endl;
	cout << "Param = " << Param << ";   KrylovSubSpace = " << Krylov << endl;

	if (qn_S > 4 || qn_S < 0)	return 1;

	cout << "\nEnter m and N: ";
	cin >> m;
	cin >> N;
	
	CInfDMRGKondoQN	DMRG(m, N, t, alpha, UL, UA, J, Krylov);
	Vector<double>	TargetState(1000, 0.0);

	double GSE = DMRG.Iterations(L, TargetState, qn_LA, qn_d, qn_S, Param);

	cout << "\nGround State Energy = " << GSE;
	cout << "; GSE / L = " << (GSE / L) << endl;

	return 0;
}
