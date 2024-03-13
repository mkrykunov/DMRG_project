//
// File InfDMRGKondoQN.h
//
// Interface of class CInfDMRGKondoQN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Kondo model
//
class CInfDMRGKondoQN
{
private:
	int		m_Size;
	int		m_m;
	double	m_t;
	double	m_alpha;
	double	m_UL;
	double	m_UA;
	double	m_J;
	double	m_EPS_DMRG;
	double	m_EPS_Diag;
	int		m_Iterations;
	int		m_Krylov;

	Matrix<double>* A1a;
	Matrix<double>* L4a;
	Matrix<double>* A1b;
	Matrix<double>* L4b;

	Matrix<double>* H1;
	Matrix<double>* H4;

	Matrix<int>*	QN_LA_1;
	Vector<double>* QN_d_1;
	Matrix<int>*	QN_LA_4;
	Vector<double>* QN_d_4;

public:
	CInfDMRGKondoQN(int, int, double, double, double, double, double, int);
	~CInfDMRGKondoQN();

	int		GetMaxSize() { return m_Size; }
	double	Iterations(int&, Vector<double>&, int, int, double, int);

private:
	int		IsExit(int, double, int);
	void	SetupStandardBlocks(void);
	int		GetSubSize(int*, int, int, double);
	void	FormSuperBlock(int*, CSparseSymmetric<double, long>&);
	double	Diagonalize_HSB(int, CSparseSymmetric<double, long>&, Vector<double>&);
	
	void	FormDensityMatrix12(int*, Vector<double>&, Vector<double>&, Matrix<int>&, Vector<double>&);
	
	void	FormDensityMatrix34(int*, Vector<double>&, Vector<double>&, Matrix<int>&, Vector<double>&);

	double	Diagonalize_ro(int, Vector<double>&, Vector<double>&, Vector<double>&,
				Matrix<int>&, Vector<double>&, Vector<long>&);

	void	Form12BlockSystem(int*, CSparseSymmetric<double, long>&, CSparseUnitary<double, long>&,
				CSparseUnitary<double, long>&, Vector<long>&);

	void	Form34BlockSystem(int*, CSparseSymmetric<double, long>&, CSparseUnitary<double, long>&,
				CSparseUnitary<double, long>&, Vector<long>&);

	void	Arrange_ro(int, Matrix<int>&, Vector<double>&, Matrix<int>&, Vector<double>&,
				Vector<double>&, Vector<long>&);
	
	void	Arrange_QN(int, int, Vector<double>&, Matrix<int>&, Vector<double>&, Matrix<int>&, Vector<double>&);

	void	Truncation(int, int, Vector<double>&, CSparse<double, long>&, Vector<double>&, Matrix<double>&);
};
