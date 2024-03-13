//
// File InfDMRG_S12_QN.h
//
// Interface of class CInfDMRG_S12_QN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Heisenberg S=1/2 model
//
class CInfDMRG_S12_QN
{
private:
	int	m_Size;
	int	m_m;
	double	m_J1;
	double	m_J2;
	double	m_Accuracy;
	double	m_EPS_DMRG;
	double	m_EPS_Diag;
	int	m_Iterations;
	
	Matrix<double>* S_up_1;
	Matrix<double>* S_up_4;
	Matrix<double>* S_z_1;
	Matrix<double>* S_z_4;

	Matrix<double>* H1;
	Matrix<double>* H4;

	Vector<double>* QN1;
	Vector<double>* QN4;

public:
	CInfDMRG_S12_QN(int, double, double, double);
	~CInfDMRG_S12_QN();

	int	GetMaxSize() { return m_Size; }
	double	Iterations(int&, Vector<double>&, double, int);

private:
	int	IsExit(int, double, int);
	void	Initialization(void);
	int	GetSubSize(int*, double);
	void	FormSuperBlock(int*, CSparseSymmetric<double, long>&, int);
	double	Diagonalize_HSB(int, CSparseSymmetric<double, long>&, Vector<double>&);
	
	void	FormDensityMatrix12(int*, Vector<double>&, Vector<double>&, Vector<double>&);

	void	FormDensityMatrix34(int*, Vector<double>&, Vector<double>&, Vector<double>&);

	double	Diagonalize_ro(int, Vector<double>&, Vector<double>&, Vector<double>&,
				Vector<double>&, Vector<long>&);

	void	Form12BlockSystem(int*, CSparseSymmetric<double, long>&,
				CSparseUnitary<double, long>&,
				CSparseUnitary<double, long>&,
				Vector<long>&, int);

	void	Form34BlockSystem(int*, CSparseSymmetric<double, long>&, CSparseUnitary<double, long>&,
				CSparseUnitary<double, long>&, Vector<long>&, int);
	
	void	Arrange_QN(int, int, Vector<double>&, Vector<double>&, Vector<double>&);

	void	Arrange_ro(int, Vector<double>&, Vector<double>&, Vector<double>&, Vector<long>&);

	void	Truncation(int, int, Vector<double>&, CSparse<double, long>&, Vector<double>&,
				Matrix<double>&);
};
