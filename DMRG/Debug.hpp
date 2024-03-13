template <class TYPE>
void Debug(int Str, int Col, Matrix<TYPE>& M, TYPE threshold)
{
	for (int i = 0; i < Str; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			if (fabs(M[i][j]) > fabs(threshold))
			{
				cout << M[i][j] << " ";
			}
			else
			{
				cout << 0 << " ";
			}
		}
		cout << endl;
	}
}

template <class TYPE>
void Debug(int Size, Vector<TYPE>& V, TYPE threshold)
{
	for (int i = 0; i < Size; i++)
	{
		if (fabs(V[i]) > fabs(threshold))
		{
			cout << V[i] << " ";
		}
		else
		{
			cout << 0 << " ";
		}
	}
	cout << endl;
}

template <class TYPE>
void Debug(int Str, int Col, Vector<TYPE>& M, TYPE threshold)
{
	for (int i = 0; i < Str; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			if (fabs(M[i * Col + j]) > fabs(threshold))
			{
				cout << M[i * Col + j] << " ";
			}
			else
			{
				cout << 0 << " ";
			}
		}
		cout << endl;
	}
}

template <class TYPE>
void Debug(int Size, Vector<TYPE>& V)
{
	for (int i = 0; i < Size; i++)
	{
		cout << V[i] << " ";
	}
	cout << endl;
}

template <class TYPE>
void Debug(int Size, Vector<TYPE>& V1, Vector<TYPE>& V2)
{
	for (int i = 0; i < Size; i++)
	{
		cout << "(" << V1[i] << "," << V2[i] << ") ";
	}
	cout << endl;
}
