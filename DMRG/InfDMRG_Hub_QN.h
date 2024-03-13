//
// File InfDMRG_Hub_QN.h
//
// Interface of class CInfDMRG_Hub_QN
// for infinite system density-matrix algorthm
// (with quantum numbers)
// for the 1D Hubbard model
//
class CInfDMRG_Hub_QN
{
private:
	int		m_Size;
	int		m_m;
	double	m_t1;
	double	m_t2;
	double	m_Accuracy;
	double	m_EPS_DMRG;
	double	m_EPS_Diag;

	Matrix<double>* c1a;
	Matrix<double>* c4a;
	Matrix<double>* c1b;
	Matrix<double>* c4b;

	Matrix<double>* H1;
	Matrix<double>* H4;

	Vector<int>* QN_up_1;
	Vector<int>* QN_down_1;
	Vector<int>* QN_up_4;
	Vector<int>* QN_down_4;

public:
	CInfDMRG_Hub_QN(int, double, double, double);
	~CInfDMRG_Hub_QN();

	int GetMaxSize() { return m_Size; }
	double Iterations(int&, Vector<double>&, int, int);

private:
	void Initialization(double);
	int	 GetSubSize(int*, int, int);
	void FormSuperBlock(int*, CSparseSymmetric<double, long>&, int);
	double Diagonalize_HSB(int, CSparseSymmetric<double, long>&, Vector<double>&);
	
	void FormDensityMatrix12(int*, Vector<double>&, Vector<double>&,
							 Vector<int>&, Vector<int>&);

	void FormDensityMatrix34(int*, Vector<double>&, Vector<double>&,
							 Vector<int>&, Vector<int>&);

	double Diagonalize_ro(int, Vector<double>&, Vector<double>&, Vector<double>&,
						  Vector<int>&, Vector<int>&,  Vector<long>&);

	void Form12BlockSystem(int*, CSparseSymmetric<double, long>&,
						   CSparseUnitary<double, long>&,
						   CSparseUnitary<double, long>&, Vector<long>&, int);
	
	void Form34BlockSystem(int*, CSparseSymmetric<double, long>&,
						   CSparseUnitary<double, long>&,
						   CSparseUnitary<double, long>&, Vector<long>&, int);
	
	void Arrange_ro(int, Vector<int>&, Vector<int>&, Vector<int>&, Vector<int>&,
					Vector<double>&, Vector<long>&);
	
	void Arrange_QN(int, int, Vector<double>&, Vector<int>&, Vector<int>&,
					Vector<int>&, Vector<int>&);

	void Truncation(int, int, Vector<double>&, CSparse<double, long>&,
					Vector<double>&, Matrix<double>&);
};
