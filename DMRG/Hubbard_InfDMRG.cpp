// g++ -o Hubbard_InfDMRG.exe Hubbard_InfDMRG.cpp External.cpp InfDMRG_Hub_QN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_2_c.c Davidson.cpp
//
#include <iostream>
#include <vector>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "CSparse.hpp"
#include "InfDMRG_Hub_QN.h"


int main(void)
{
	int L = 32;
	int m = 20;

	double t1 = -1.0;
	double t2 = -1.0;
	double U = 10;

	int qn_up = 2;
	int qn_down = 2;

	FILE* stream;

	if ((stream = fopen("DAT\\Hubbard.dat", "r")) != NULL)
	{
		fscanf(stream, "%lf %lf %lf %d %d", &t1, &t2, &U, &qn_up, &qn_down);
		fclose(stream);
	}

	cout << "t1 = " << t1 << ";   t2 = " << t2 << ";   U = " << U;
	cout << ";   qn_up = " << qn_up << ";   qn_down = " << qn_down << endl;

	cout << "\nEnter m: ";
	cin >> m;

	CInfDMRG_Hub_QN DMRG(m, t1, t2, U);
	
	Vector<double> TargetState;

	double GroundStateEnergy = DMRG.Iterations(L, TargetState, qn_up, qn_down);

	cout << "\nGround State Energy = " << GroundStateEnergy;
	cout << "; GSE / L = " << (GroundStateEnergy / L) << endl;

	return 0;
}
