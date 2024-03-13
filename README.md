# DMRG_project
This is one of my PhD projects. It contains the code for three models: Heisenberg, Hubbard and Kondo.

To compile the infinite DMRG algorithm for the Heisenberg model:

g++ -o Heisenberg_InfDMRG.exe Heisenberg_InfDMRG.cpp External.cpp InfDMRG_S12_QN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_c.c Davidson.cpp

To run Heisenberg_InfDMRG.exe reads the parameters from Dat/Heisenberg.dat
When it starts, it askes for m, enter for example 32
It should converge after 92 iterations.

To compile the infinite DMRG algorithm for the Hubbard model:

g++ -o Hubbard_InfDMRG.exe Hubbard_InfDMRG.cpp External.cpp InfDMRG_Hub_QN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_2_c.c Davidson.cpp

To run Hubbard_InfDMRG.exe reads the parameters from Dat/Hubbard.dat
When it starts, it askes for m, enter for example 64
It should converge after 92 iterations.

To compile the infinite DMRG algorithm for the Kondo model:

g++ -o Kondo_InfDMRG.exe Kondo_InfDMRG.cpp External.cpp InfDMRGKondoQN.cpp tred2_c.c tqli_c.c sort_c.c sort_qn_3_c.c Davidson.cpp

To run Kondo_InfDMRG.exe reads the parameters from Dat/Kondo.dat
When it starts, it askes for m and N, enter for example 64 4
It should converge after 100 iterations.
