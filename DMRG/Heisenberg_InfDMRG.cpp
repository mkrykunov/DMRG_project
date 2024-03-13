// g++ -o Heisenberg_InfDMRG.exe Heisenberg_InfDMRG.cpp External.cpp InfDMRG_S12_QN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_c.c Davidson.cpp
// 
#include <iostream>
#include <vector>

#include "External.h"
#include "Linalg.hpp"
#include "EigenSys.hpp"
#include "CSparse.hpp"

#include "Spin2.h"
//#include "Spin3.h"

#include "InfDMRG_S12_QN.h"


int main(void)
{
	int L = 32;
	int m = 20;
	int Param = 0;
	
	double J1 = 1.0;
	double J2 = 1.0;

#if __SPIN == 2
	double S = 0.5;
	double qn_Sz = 0.5;
#elif __SPIN == 3
	double S = 1.0;
	double qn_Sz = 1.0;
#endif
	
//	double S = 0.5;
//	double qn_Sz = 0.5;

	FILE* stream;

	if ((stream = fopen("DAT\\Heisenberg.dat", "r")) != NULL)
	{
		fscanf(stream, "%lf %lf %lf %d", &J1, &J2, &qn_Sz, &Param);
		fclose(stream);
	}

	std::cout << "J1 = " << J1 << ";   J2 = " << J2 << ";   qn_Sz = " << qn_Sz << std::endl;
	std::cout << "Param = " << Param << std::endl;

	std::cout << "\nEnter m: ";
	std::cin >> m;

	CInfDMRG_S12_QN DMRG(m, J1, J2, S);
	
//	int Size = DMRG.GetMaxSize();

	Vector<double> TargetState(1000, 0.0);

	double GroundStateEnergy = DMRG.Iterations(L, TargetState, qn_Sz, Param);

	std::cout << "\nGround State Energy = " << GroundStateEnergy;
	std::cout << "; GSE / L = " << (GroundStateEnergy / L) << std::endl;

	return 0;
}
