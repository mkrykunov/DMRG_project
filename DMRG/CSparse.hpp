// Sparse storage for matrix
//
template <class TYPE1, class TYPE2>
class CSparse
{
public:
	vector<TYPE1, allocator<TYPE1> > sa;

	virtual double GetSparsity(void) { return 0.0; }
	virtual void Multiply(TYPE2 N1, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b) {}
	virtual void Multiply_TR(TYPE2 N1, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b) {}
};

// Sparse storage for symmetric matrix
//
template <class TYPE1, class TYPE2>
class CSparseSymmetric : public CSparse<TYPE1, TYPE2>
{
public:
	vector<TYPE2, allocator<TYPE2> > ir;
	vector<TYPE2, allocator<TYPE2> > ic;

public:
	void resize(int i1, int i2, int i3)
	{
		this->sa.resize(i1, (TYPE1)0);
		ir.resize(i2, (TYPE2)0);
		ic.resize(i3, (TYPE2)0);
	}

	void reserve(int i1, int i2, int i3)
	{
		this->sa.reserve(i1);
		ir.reserve(i2);
		ic.reserve(i3);
	}

	inline void push_back(TYPE1 a1, TYPE2 a2, TYPE2 a3)
	{
		if (a2 > a3)
		{
			return;
		}
		else if (a2 == a3)
		{
			this->sa[a2] = a1;
		}
		else if (a1 != (TYPE1)0)
		{
			this->sa.push_back(a1);
			ir.push_back(a2);
			ic.push_back(a3);
		}
	}

	double GetSparsity(void)
	{
		if (ir.size() != ic.size())
		{
			return 0.0;
		}
		else if (this->sa.size() != 0 && ir.size() != 0)
		{
			TYPE2 Size = this->sa.size() - ir.size();

			return (Size < 0) ? 0.0 : 100.0 * this->sa.size() / (0.5 * Size * (Size + 1));
		}
		else
		{
			return 0.0;
		}
	}

	void Multiply(TYPE2 N1, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b)
	{
		if (ir.size() != ic.size())	return;

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < M; j++)
			{
				c[i * M + j] = this->sa[i] * b[i * N1 + j];
			}
		}

		for (int ii = 0; ii < ir.size(); ii++)
		{
			int i = ir[ii];
			int k = ic[ii];

			double aik = this->sa[N1 + ii];

			for (int j = 0; j < M; j++)
			{
				c[i * M + j] += aik * b[k * N1 + j];
				c[k * M + j] += aik * b[i * N1 + j];
			}
		}
	}

	void Multiply_TR(TYPE2 N, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b)
	{
		if (ir.size() != ic.size())	return;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				c[i * M + j] = this->sa[i] * b[j * N + i];
			}
		}

		for (int ii = 0; ii < ir.size(); ii++)
		{
			int i = ir[ii];
			int k = ic[ii];

			double aik = this->sa[N + ii];

			for (int j = 0; j < M; j++)
			{
				c[i * M + j] += aik * b[j * N + k];
				c[k * M + j] += aik * b[j * N + i];
			}
		}
	}
};

// Sparse storage for unitary matrix
//
template <class TYPE1, class TYPE2>
class CSparseUnitary : public CSparse<TYPE1, TYPE2>
{
public:
	vector<TYPE2, allocator<TYPE2> > ija;

public:
	void resize(int i1, int i2)
	{
		this->sa.resize(i1, (TYPE1)0);
		ija.resize(i2, (TYPE2)0);
	}

	void reserve(int i1, int i2)
	{
		this->sa.reserve(i1);
		ija.reserve(i2);
	}

	void push_back(TYPE1 a1, TYPE2 a2, TYPE2 a3)
	{
		if (a2 == a3)
		{
			this->sa[a2] = a1;
		}
		else if (a1 != (TYPE1)0)
		{
			this->sa.push_back(a1);
			ija.push_back(a3);
		}
	}

	double GetSparsity(void)
	{
		if (this->sa.size() != 0 && ija.size() != 0)
		{
			TYPE2 Size = ija[0] - 1;

			return (Size < 0) ? 0.0 : 100.0 * this->sa.size() / (Size * Size);
		}
		else
		{
			return 0.0;
		}
	}

	void Multiply(TYPE2 N, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b)
	{
		for (TYPE2 i = 0; i < N; i++)
		{
			for (TYPE2 j = 0; j < M; j++)
			{
				c[i * M + j] = this->sa[i] * b[i * N + j];

				for (TYPE2 k = ija[i]; k < ija[i+1]; k++)
				{
					c[i * M + j] += this->sa[k] * b[ija[k] * N + j];
				}
			}
		}
	}

	void Multiply_TR(TYPE2 N, TYPE2 M, Vector<TYPE1>& c, Vector<TYPE1>& b)
	{
		for (TYPE2 i = 0; i < N; i++)
		{
			for (TYPE2 j = 0; j < M; j++)
			{
				c[i * M + j] = this->sa[i] * b[j * N + i];

				for (TYPE2 k = ija[i]; k < ija[i+1]; k++)
				{
					c[i * M + j] += this->sa[k] * b[j * N + ija[k]];
				}
			}
		}
	}
};
