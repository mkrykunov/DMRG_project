/////////////////////////////////////////////////////////////////////////////
// Interfaces of Linalg classes

using namespace std;

const unsigned short CONSTANT  = 1;
const unsigned short TEMPORAL  = 2;
const unsigned short SYMMETRIC = 4;
const unsigned short FULLDIAG  = 8;
const unsigned short EXTERNAL  = 16;

const char mem_fail[] = "Memory allocation failure!";

template <class TYPE> class Matrix;
template <class TYPE> class SparseMatrix;

/////////////////////////////////////////////////////////////////////////////
// Interface of the Vector class

#ifndef  __VECTOR_H_
#define  __VECTOR_H_

template <class TYPE> class Vector 
{
// Attributes
protected:
	unsigned short m_Flag;
	long	Size;
	TYPE*   Vec; 

// Construction / Destruction
public:
	Vector();
	Vector(Vector<TYPE>&);
	Vector(long, unsigned short Stat = CONSTANT);
	Vector(long, TYPE, unsigned short Stat = CONSTANT);
	void Realloc(long);
	~Vector();

// Operations
public: 
	void	operator=(Vector<TYPE>*);  
	Vector<TYPE>&	operator=(Vector<TYPE>&);

	TYPE&	operator[](long pos) { return  Vec[pos]; }
	TYPE&	operator()(long pos) { return  Vec[pos]; }

	void Read( const char* );
	void Write( const char*, int );
	template <class U> friend ostream& operator << (ostream&, Vector<U>&);
    
	unsigned short getFlag() const { return  m_Flag; }
	long	getSize() const { return  Size; }
	TYPE*	getVec() const { return  Vec; }
	TYPE	Norma();
	void	AddMult(SparseMatrix<TYPE>&, Vector<TYPE>&, long);	// V += M * V;
    
	template<class U> friend Vector<U>& operator*(Matrix<U>&, Vector<U>&);
	template<class U> friend Vector<U>& operator+(Vector<U>&, Vector<U>&);
	template<class U> friend Vector<U>& operator-(Vector<U>&, Vector<U>&);
	template<class U> friend Vector<U>& operator*(U&, Vector<U>&);   
	template<class U> friend U          operator*(Vector<U>&, Vector<U>&); 
};


template <class TYPE> Vector<TYPE>::Vector()
{
	m_Flag = CONSTANT;
	Size   = 0; 
	Vec    = NULL;
}   

template <class TYPE> Vector<TYPE>::Vector(Vector<TYPE>& aVector)
{
	m_Flag = aVector.m_Flag;
	Size = aVector.Size; 

	if (!(Vec = new TYPE [Size]))
		Error(mem_fail);
    
	for (long i = 0; i < Size; i++)
		Vec[i] = (TYPE)aVector[i];  
}   

template <class TYPE> Vector<TYPE>::Vector(long arg, unsigned short Stat)
{
	m_Flag = Stat;
	Size = arg; 

	if (!(Vec = new TYPE [Size]))
		Error(mem_fail);
    
	for (long i = 0; i < Size; i++)	Vec[i] = (TYPE)0.0;  
}   

template <class TYPE> Vector<TYPE>::Vector(long arg, TYPE num, unsigned short Stat)
{
	m_Flag = Stat;
	Size = arg;

	if (!(Vec = new TYPE [Size]))
		Error(mem_fail);
    
	for (long i = 0; i < Size; i++)	Vec[i] = num;
}   

template <class TYPE> void Vector<TYPE>::Realloc(long arg)
{
	if (Size != 0)	delete [] Vec;

	Size = arg; 

	if (!(Vec = new TYPE [Size]))
		Error(mem_fail);
    
	for (long i = 0; i < Size; i++)	Vec[i] = (TYPE)0.0;  
}   

template <class TYPE> Vector<TYPE>::~Vector()    
{
	if ((m_Flag & EXTERNAL) != EXTERNAL)
	{
		delete [] Vec;
	}
}  

#endif // __VECTOR_H_

/////////////////////////////////////////////////////////////////////////////
// Interface of the Matrix class

#ifndef  __MATRIX_H_
#define  __MATRIX_H_ 

template <class TYPE> class Matrix
{
// Attributes
private:
	unsigned short m_Flag;
	long	m_nStrings;
	long	m_nColumns;
	TYPE**  Mat;

// Construction / Destruction
public:             
	Matrix(long, long, unsigned short Stat = CONSTANT);
	Matrix(Matrix<TYPE>&);
	Matrix(long, long, TYPE, unsigned short Stat = CONSTANT);
	Matrix(SparseMatrix<TYPE>&, unsigned short /*Flag*/);
	void Realloc(long, long);
	~Matrix();

// Operations
public:       
	void operator=(Matrix<TYPE>*);
	Matrix<TYPE>& operator=(Matrix<TYPE>&);

	TYPE*   operator[](long pos) { return  Mat[pos]; }
	TYPE&   operator()(long i, long j) { return  Mat[i][j]; }
	
	void Read(const char*); 
	void Write(const char*, int); 
	void Write(int);
	template<class U> friend ostream& operator << (ostream&, Matrix<U>&);

	long	getStrings() const { return  m_nStrings; }
	long	getColumns() const { return  m_nColumns; }
	TYPE**	getMat() const { return  Mat; }
    
	template<class U> friend U Spur(Matrix<U>&);
    
	template<class U> friend Matrix<U>& operator+(Matrix<U>&, Matrix<U>&); 
	template<class U> friend Matrix<U>& operator-(Matrix<U>&, Matrix<U>&);
	template<class U> friend Matrix<U>& operator*(Matrix<U>&, Matrix<U>&);
	template<class U> friend Vector<U>& operator*(Matrix<U>&, Vector<U>&);
	template<class U> friend Matrix<U>& operator*(U&, Matrix<U>&);   
};

template <class TYPE> Matrix<TYPE>::Matrix(long arg1, long arg2, unsigned short Stat)
{
	m_Flag = Stat;
	m_nStrings = arg1;
	m_nColumns = arg2;

	if (!(Mat = new TYPE*[m_nStrings]))
		Error(mem_fail); 
       
	for (long i = 0; i < m_nStrings; i++)  
		if (!(Mat[i] = new TYPE[m_nColumns]))
			Error(mem_fail);

 	for (long i = 0; i < m_nStrings; i++)  
		for (long j = 0; j < m_nColumns; j++)   Mat[i][j] = 0.0;
}

template <class TYPE> Matrix<TYPE>::Matrix(Matrix<TYPE>& aMatrix)
{
	m_Flag = aMatrix.m_Flag;
	m_nStrings = aMatrix.m_nStrings;
	m_nColumns = aMatrix.m_nColumns;

	if (!(Mat = new TYPE* [m_nStrings]))
		Error(mem_fail); 
       
	for (long i = 0; i < m_nStrings; i++)  
		if (!(Mat[i] = new TYPE [m_nColumns]))
			Error(mem_fail);

	for (long i = 0; i < m_nStrings; i++)  
		for (long j = 0; j < m_nColumns; j++)   Mat[i][j] = aMatrix[i][j];
}

template <class TYPE> Matrix<TYPE>::Matrix(long arg1, long arg2, TYPE num, unsigned short Stat)
{
	m_Flag = Stat;
	m_nStrings = arg1;
	m_nColumns = arg2;

	if (!(Mat = new TYPE*[m_nStrings]))
		Error(mem_fail); 
       
	for (long i = 0; i < m_nStrings; i++)  
		if (!(Mat[i] = new TYPE[m_nColumns]))
			Error(mem_fail);

	for (long i = 0; i < m_nStrings; i++)  
		for (long j = 0; j < m_nColumns; j++)   Mat[i][j] = num;
}

template <class TYPE> Matrix<TYPE>::Matrix(SparseMatrix<TYPE>& arg, unsigned short Stat)
{
	m_Flag = CONSTANT;
	m_nStrings = arg.getSize();
	m_nColumns = arg.getSize();    
    
	if (!(Mat = new TYPE*[m_nStrings]))
	       Error(mem_fail); 
       
	for (long i = 0; i < m_nStrings; i++)
		if (!(Mat[i] = new TYPE[m_nColumns]))
			Error(mem_fail);

	for (long j = 0; j < m_nColumns; j++)  
		for (long i = 0; i < m_nStrings; i++)
			Mat[i][j] = 0.0;

	for (long k = 0; k < arg.getNEL(); k++)
	{
		long i = arg.r(k);
		long j = arg.c(k);
		Mat[i][j] = arg[k];
		if (Stat == SYMMETRIC)
			Mat[j][i] = Mat[i][j];
	}
}

template <class TYPE> void Matrix<TYPE>::Realloc(long arg1, long arg2)
{
	if (m_nColumns != 0)
	{
		for (long i = 0; i < m_nStrings; i++)
			delete [] Mat[i];
	}

	if (m_nStrings != 0)
		delete  Mat;

	m_nStrings = arg1;
	m_nColumns = arg2;

	if (!(Mat = new TYPE*[m_nStrings]))
		Error(mem_fail); 
       
	for (long i = 0; i < m_nStrings; i++)  
		if (!(Mat[i] = new TYPE[m_nColumns]))
			Error(mem_fail);

	for (long i = 0; i < m_nStrings; i++)  
		for (long j = 0; j < m_nColumns; j++)   Mat[i][j] = (TYPE)0.0;
}   

template <class TYPE> Matrix<TYPE>::~Matrix()
{
	for (long i = 0; i < m_nStrings; i++)
		delete [] Mat[i];

	delete  Mat;
}

#endif // __MATRIX_H_

/////////////////////////////////////////////////////////////////////////////
// Sparse.h : interface of the SparseMatrix class

#ifndef  __SPARSE_H_
#define  __SPARSE_H_

template <class TYPE> class SparseMatrix : public Vector<TYPE>
{
// Attributes
private:
	long  NEL;
	long* IR;
	long* IC;

// Construction / Destruction
public:
	SparseMatrix(long, long, const char Status = 0);
	SparseMatrix(Matrix<TYPE>*, const char);
	SparseMatrix(long, long, TYPE*, long*, long*, const char);
	~SparseMatrix();
    
// Operations
public:
	long&	r(long arg) { return IR[arg]; }
	long&	c(long arg) { return IC[arg]; }
	TYPE&	push_back(long, long, long);
	TYPE&	add(long, long, long);
	void	add(long, long, TYPE, long);
	TYPE	e(long, long);
	long	getNEL() const { return  NEL; }
	long*	getIR() const { return  IR; }
	long*	getIC() const { return  IC; }

//	void Read( const char* ); 
	void Write(const char*, int);

	long DiagToBegin();
	void ConvertToFullDiag();
	void __ConvertToFullDiag();
	void ConvertFromFullDiag();
};

template <class TYPE>
SparseMatrix<TYPE>::SparseMatrix(long arg1, long arg2, const char Status)
				 : Vector<TYPE>(arg1 + arg2)
{
	this->m_Flag = Status;
	this->Size = arg1;
	NEL = ((Status & FULLDIAG) == 0) ? arg2	: 0;
    
	if (!(IR = new long [arg2]))	Error(mem_fail);
	if (!(IC = new long [arg2]))	Error(mem_fail);
}

template <class TYPE>
SparseMatrix<TYPE>::SparseMatrix(Matrix<TYPE>* arg, const char Status)
{
	if ((Status & SYMMETRIC) == 0)
		Error("Unknown Status");

	this->m_Flag = Status;
	this->Size = arg->getStrings();
	NEL = 0;
    
	for (long i = 0; i < this->Size; i++)
		for (long j = i; j < this->Size; j++)
		    if ((*arg)[i][j] != 0.0)    NEL++;
		    
	if (!(this->Vec = new TYPE [NEL]))	Error(mem_fail);
	if (!(IR = new long [NEL]))		Error(mem_fail);
	if (!(IC = new long [NEL]))		Error(mem_fail);
    
	long Counter = 0;

	for (long i = 0; i < this->Size; i++)
	{
		for (long j = i; j < this->Size; j++)
		{
			if ((*arg)[i][j] != 0.0)
			{
				this->Vec[Counter] = (*arg)[i][j];
				IR[Counter] = i;
				IC[Counter] = j;
				Counter++;
			}
		}
	}
}

template <class TYPE>
SparseMatrix<TYPE>::SparseMatrix(long aSize, long aNEL, TYPE* sa, long* ir,
								 long* ic, const char Status)
{
	if ((Status & SYMMETRIC) == 0)
		Error("Unknown Status");

	this->m_Flag = Status | EXTERNAL;
	this->Size = aSize;
	NEL = aNEL;

	this->Vec = sa;
	IR = ir;
	IC = ic;
}

template <class TYPE>
SparseMatrix<TYPE>::~SparseMatrix()
{
	if ((this->m_Flag & EXTERNAL) != EXTERNAL)
	{
		delete [] IR;
		delete [] IC;
	}
}

template <class TYPE>
TYPE& SparseMatrix<TYPE>::push_back(long Row, long Column, long pos)
{
	NEL++;

	IR[pos] = Row;
	IC[pos] = Column;

	if ((this->m_Flag & FULLDIAG) != 0)
		pos += this->Size;

	return  this->Vec[pos];
}


#endif // __SPARSE_H_

/////////////////////////////////////////////////////////////////////////////
// Globals

template <class TYPE> TYPE Norming(long, TYPE*);

template <class TYPE> TYPE** allocate(TYPE** ptr, int N, int M)
{
	TYPE** pointer = new TYPE* [N];

	if (!pointer)	Error(mem_fail);

	for (int i = 0; i < N; i++)
	{
		if (!(pointer[i] = new TYPE [M]))	Error(mem_fail);
	}

	return pointer;
}


template <class TYPE> void deallocate(TYPE** pointer, int N)
{
	for (int i = 0; i < N; i++)
		delete [] pointer[i];

	if (pointer)	delete [] pointer;
}

/////////////////////////////////////////////////////////////////////////////
// Matrix Globals

template <class TYPE> TYPE Scalar(Matrix<TYPE>&, Vector<TYPE>&);	// <xHx> - function
template <class TYPE> TYPE Det( Matrix<TYPE>&);	// Determinant of Matrix
template <class TYPE> Matrix<TYPE>& Inverse(Matrix<TYPE>&);	// Gausse Inversion of Matrix

/////////////////////////////////////////////////////////////////////////////
// Vector Globals

template <class TYPE> TYPE MaxAbsElement(Vector<TYPE>&, Vector<TYPE>&);	// Maximal Value of |a(i)-b(i)|
template <class TYPE> TYPE MaxAbsElement(Matrix<TYPE>&, Matrix<TYPE>&);	//          --//--
template <class TYPE> TYPE MaxAbsElement(long, TYPE*, TYPE*);	        //          --//--
template <class TYPE> TYPE Norming(Vector<TYPE>&);	// Norming of Vector

/////////////////////////////////////////////////////////////////////////////
// SparseMatrix Globals

template <class TYPE> TYPE Scalar(SparseMatrix<TYPE>&, TYPE*, long);			// <xHx> - function
template <class TYPE> void AddMult(TYPE*, SparseMatrix<TYPE>&, TYPE*, long);	// V += M * V;

