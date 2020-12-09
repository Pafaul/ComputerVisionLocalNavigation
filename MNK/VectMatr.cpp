#include "VectMatr.h"

TVector::TVector(int comesize) :size(comesize), value(nullptr)
{
	this->value = new double[comesize];
};
// ����������� �����������
TVector::TVector(const TVector& rvalue) : size(0), value(nullptr) 
{
	this->size = rvalue.size;
	value = new double[rvalue.size];
	for (int i = 0; i < rvalue.size; i++)
	{
		this->value[i] = rvalue.value[i];
	}
}
TVector::~TVector()
{
	if (value) 
	{
		delete[] value;
		size = 0;
		value = nullptr;
	}
};
void TVector::operator = (const TVector& comevect)
{
	//�������� ������ ������� � ���������� ����� ���������������
	if (this->value)
	{
		delete[] this->value;
		this->size = 0;
		this->value = nullptr;
	}
	//��� �������������� ������� ����������
	if (this->size != comevect.size)
	{
		this->resize(comevect.size);
	}
	//������ ������������
	for (int i = 0; i < size; i++)
	{
		this->value[i] = comevect[i];
	}
};
TVector TVector::operator + (const TVector& comevect) const
{
	TVector sub(this->size);
	//���� ���������� �� �������� ������ ������� ��� ����� ��������
	for (int i = 0; i < this->size; i++)
	{
		sub[i] = this->value[i] + comevect[i];
	}

	return sub;
};
TVector TVector::operator - (const TVector& comevect) const
{
	TVector sub(this->size);
	//���� ���������� �� �������� ������ ������� ��� ����� ��������
	for (int i = 0; i < this->size; i++)
	{
		sub[i] = this->value[i] - comevect[i];
	}

	return sub;
}
double TVector::operator * (const TVector& comevect) const
{
	double sub = 0;
	//���� ����������
	for (int i = 0; i < size; i++)
	{
		sub += this->value[i] * comevect.value[i];
	}
	return sub;
};
TVector TVector::operator * (const double& comeval) const
{
	TVector sub(size);
	for (int i = 0; i < size; i++)
	{
		sub[i] = this->value[i] * comeval;
	}
	return sub;
};
TVector TVector::operator ^ (const TVector& comevect) const
{
	TVector sub(this->size);
	//���������� ���������� ������������
	sub[0] = this->value[1] * comevect[2] - this->value[2] * comevect[1];
	sub[1] = this->value[2] * comevect[0] - this->value[0] * comevect[2];
	sub[2] = this->value[0] * comevect[1] - this->value[1] * comevect[0];

	return sub;
};
void TVector::add(double comeval[])
{ 
	for (int i = 0; i < size; i++)
	{
		this->value[i] = comeval[i];
	};
}
void TVector::add(std::vector<double> vec)
{
	for (int i = 0; i < size; i++)
	{
		this->value[i] = vec[i];
	};
}
;
double TVector::mod() const
{
	double sum = 0;
	for (int i = 0; i < size; i++)
	{
		sum += std::pow(value[i], 2);
	}
	return sqrt(sum);
};
TVector TVector::norm() const
{
	TVector sub(size);
	double inleng = mod();
	for (int i = 0; i < size; i++)
	{
		sub.value[i] = value[i] / inleng;
	}
	return sub;
}
//������� ��������� �����
TVector TVector::negative()
{
	TVector sub(size);
	for (int i = 0; i < size; i++)
	{
		sub.value[i] = -value[i];
	}
	return sub;
}
//����������
void TVector::resize(int new_size)
{
	TVector sub(new_size);
	int min_size = new_size < this->size ? new_size : this->size;
	if (this->value)
	{
		for (int i = 0; i < min_size; i++) sub.value[i] = this->value[i];
		delete[] value;
		this->size = 0;
		this->value = nullptr;
	}
	this->value = new double[new_size];
	this->size = new_size;
	if (sub.value) for (int i = 0; i < min_size; i++) this->value[i] = sub[i];
}
// �������� - ������� �����
TVector TVector::operator - () const
{
	TVector sub(size);
	for (int i = 0; i < size; i++)
	{
		sub.value[i] = -value[i];
	}
	return sub;
};
//������� ������� ������ �������� ��� �� �������� ���� ��� ������ ������� �������
TVector TVector::Rodrig(double phi, const TVector& axis)
{
	double phi_gr = (phi / 180) * 3.1416;
	TVector V(3);
	TVector e = axis.norm();
	V = (e * *this) * (1 - cos(phi_gr)) * e + (e ^ *this) * sin(phi_gr) + *this * cos(phi_gr);
	return V;
}
//																		TMatrix


TMatrix::TMatrix(int x, int y) :high(x), lengh(y), value(nullptr)
{
	this->value = new double*[x];
	for (int i = 0; i < x; i++)
	{
		this->value[i] = new double[y];
	}
};
//����������
TMatrix::~TMatrix()
{
	if (value)
	{
		for (int i = 0; i < high; i++)
		{
			delete[] value[i];
			value[i] = nullptr;
		}
		delete[] value;                                    
		high = 0;
		lengh = 0;
		value = nullptr;
	}
};
// ����������� �����������
TMatrix::TMatrix(const TMatrix& rvalue) :high(0), lengh(0), value(nullptr)
{
	this->high = rvalue.high;
	this->lengh = rvalue.lengh;
	value = new double*[rvalue.high];
	for (int i = 0; i < rvalue.high; i++)
	{
		this->value[i] = new double[rvalue.lengh];
	};
	for (int i = 0; i < rvalue.high; i++)
	{
		for (int j = 0; j < rvalue.lengh; j++)
		{
			this->value[i][j] = rvalue.value[i][j];
		}
	}
};
//������������
void TMatrix::operator = (const TMatrix& comematr)
{
	//�������� ������ ������� � ���������� ����� ���������������
	if (this->value)
	{
		for (int i = 0; i < high; i++)
		{
			delete[] value[i];
		}
		delete[] value;
		high = 0;
		lengh = 0;
		value = nullptr;
	}
	//��� �������������� ������ ����������
	if (this->high != comematr.high || this->lengh != comematr.lengh)
	{
		this->resize(comematr.high, comematr.lengh);
	}
	//������ ������������
	for (int i = 0; i < comematr.high; i++)
	{
		this->value[i] = new double[comematr.lengh];
		for (int j = 0; j < comematr.lengh; j++)
		{
			this->value[i][j] = comematr.value[i][j];
		}
	}
};
TMatrix TMatrix::operator + (const TMatrix &comematr)
{
	TMatrix sub(this->high, this->lengh);
	//���� ���������� �� �������� ������ ������� ��� ����� ��������
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			sub(i, j) = this->value[i][j] + comematr(i, j);
		}
	}
	
	return sub;
};
TMatrix TMatrix::operator - (const TMatrix& comematr)
{
	TMatrix sub(this->high, this->lengh);
	//���� ���������� �� �������� ������ ������� ��� ����� ��������
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			sub(i, j) = this->value[i][j] - comematr(i, j);
		}
	}
	
	return sub;
}
//������������ ������
TMatrix TMatrix::operator * (const TMatrix &comematr)//����� ��������� ��� ������������ ����������
{
	TMatrix sub(this->high, comematr.lengh);
	if (comematr.high != this->lengh)
	{
		sub(0, 0) = 0;
		return sub;
	};
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < comematr.lengh; j++)
		{
			sub(i, j) = 0;
			for (int n = 0; n < this->lengh; n++)
			{
				sub(i, j) += this->value[i][n] * comematr.value[n][j];
			}
		}
	}
	return sub;
};
TVector TMatrix::operator * (const TVector&comevect)
{
	TMatrix sub(this->high, comevect.length());
	TVector sub2(this->high);
	if (comevect.length() != this->lengh)
	{
		sub2[0] = 0;
		return sub2;
	};
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < comevect.length(); j++)
		{
			sub(i, j) = 0;
			for (int n = 0; n < this->lengh; n++)
			{
				sub(i, j) += this->value[i][n] * comevect[n];
			}
		}
	}
	for (int i = 0; i < this->high; i++)
	{
		sub2[i] = sub(i, 0);
	}
	return sub2;
}
//��������� �� �����
TMatrix TMatrix::operator * (const double& comeval)
{
	TMatrix sub(this->high, this->lengh);
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			sub(i, j) = this->value[i][j] * comeval;
		}
	}
	return sub;
};

int TMatrix::search(double what, bool match, unsigned int &uI, unsigned int &uJ, unsigned int starti, unsigned int startj) 
{
	// ����� � ������� a[m][n] �������� � ��������� ��������� what
	// ������������� ��� ����� ������ � ������� uI, uJ, ���� ������� ������.
	// match - ������ ������ ������� ��� �������� �� ����������.
	// ����� 0 - �� �������, �� 0 - �������
	if ((!lengh) || (!high)) return 0;
	if ((starti >= high) || (startj >= lengh)) return 0;
	for (unsigned int i = starti; i < high; i++)
		for (unsigned int j = startj; j < lengh; j++) {
			if (match == true) {
				if (value[i][i] == what) {
					uI = i; uJ = j; return 1;
				}
			}
			else if (value[i][j] != what) {
				uI = i; uJ = j; return 1;
			}
		}
	return 0;
}
//����� �����
void TMatrix::swaprows(unsigned int x1, unsigned int x2) {
	//������ � ������� a[n][m] ������ � �������� x1 � x2 �������
	if ((!lengh) || (!high)) return;
	if ((x1 >= lengh) || (x2 >= lengh) || (x1 == x2)) return;
	double tmp;
	for (unsigned int x = 0; x < high; x++) {
		tmp = value[x1][x];
		value[x1][x] = value[x2][x];
		value[x2][x] = tmp;
	}
	return;
};
//����� ��������
void TMatrix::swapcolumns(unsigned int x1, unsigned int x2) {
	//������ � ������� a[n][m] ������� � �������� x1 � x2 �������
	if ((!lengh) || (!high)) return;
	if ((x1 >= high) || (x2 >= high) || (x1 == x2)) return;
	double tmp;
	for (unsigned int x = 0; x < lengh; x++) {
		tmp = value[x][x1];
		value[x][x1] = value[x][x2];
		value[x][x2] = tmp;
	}
	return;
};
//�����������
double TMatrix::det()
{
	TMatrix sub(high, lengh);
	for (int i = 0; i < high; i++)
	{
		for (int j = 0; j < lengh; j++)
		{
			sub.value[i][j] = this->value[i][j];
		}
	}
	if (high = lengh)
	{
		unsigned int m = high;
		if (m == 0) return 0;
		if (m == 1) return sub.value[0][0];
		if (m == 2) return (sub.value[0][0] * sub.value[1][1] - sub.value[1][0] * sub.value[0][1]);
		bool sign = false; // ����� ����� ������������. �� ��������� - ���
		double det = 1; // ������������
		double tmp;
		unsigned int x, y;
		for (unsigned int i = 0; i < high; i++) { // ���� �� ���� ������� ���������
			if (sub.value[i][i] == 0) { // ���� ������� �� ��������� ����� 0, �� ���� ��������� ������� � �������
				if (!search(0, false, y, x, i, i)) return 0; // ���� ��� �������� �������, �� ���. = 0
				if (i != y) { // ������ i-�� ������ � y-��
					swaprows(i, y);
					sign = !sign;
				}
				if (i != x) { // ������ i-�� ������� � x-��
					swapcolumns(i, x);
					sign = !sign;
				}
				// ����� �������, � a[i][i], ������ ��������� �������.
			}
			// ������� ������� a[i][i] �� ������������
			det *= sub.value[i][i];
			tmp = sub.value[i][i];
			for (x = i; x < m; x++) {
				sub.value[i][x] = sub.value[i][x] / tmp;
			}
			// ����� ������� a[i][i] ������ ����� 1
			// �������� ��� �������� ������� ��� (i, i)-��,
			// ��� ������ ��������� � ���. �������������
			for (y = i + 1; y < high; y++) {
				tmp = sub.value[y][i];
				for (x = i; x < m; x++)
					sub.value[y][x] -= (sub.value[i][x] * tmp);
			}
		}
		if (sign) return det*(-1);
		return det;
	}
}
//�������� ����������
bool TMatrix::sylver()
{
	TMatrix sub(this->high - 1, this->lengh - 1);
	//�������� ���� �������
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			int m = 0, n = 0;
			//������ - �������
			if (j == i)
			{
				//���������� sub
				for (int i2 = 0; (i2 < this->high); i2++)
				{
					if (i2 != i)
					{
						n = 0;
						for (int j2 = 0; (j2 < this->lengh); j2++)
						{
							if (j2 != j)
							{
								sub(m, n) = value[j2][i2];
								n++;
							}
						}
						m++;
					}
				}
			}
			if (sub.det() <= 0) return false;
		}
	}
	return true;
}
// ������������� ��������
void TMatrix::add(double* comeval[])
{
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			this->value[i][j] = comeval[i][j];
		}
	}
};
void TMatrix::add(std::vector<double> vec, int length)
{
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			this->value[i][j] = vec[i*length + j];
		}
	}
}
TMatrix TMatrix::transp()
{
	TMatrix sub(this->lengh, this->high);

	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			sub.value[j][i] = this->value[i][j];
		}
	}
	return sub;
}
void TMatrix::Check_diag_zero(TMatrix& sub)
{
	double null = LDBL_EPSILON;

	for (int i = 0; i < this->high; i++)
	{
		if ((this->value[i][i] <= null) && (this->value[i][i] >= -null))
		{
			for (int k = 0; k < this->high; k++)
			{
				if ((this->value[k][i] >= null) || (this->value[k][i] <= -null))//������ �� ����� 0
				{
					if ((this->value[i][k] >= null) || (this->value[i][k] <= -null))//���������� ������ ����� �� 0
					{
						for (int j = 0; j < this->lengh; j++)
						{
							//�������� �������
							this->swaprows(i, k);
							//��������� �������
							sub.swaprows(i, k);
						}
					}
				}
			}
		}
	}
};
void TMatrix::gaus()
{ 
	//�������� ��� ������� ���������� � �������������
	if ((this->high == this->lengh) && (this->det() != 0))
	{
		TMatrix sub = E();
		//�������� ��� �� ��������� �� 0(���������)
		this->Check_diag_zero(sub);
		//�����
		//���������� � ����������������� ����
		for (int i = 1; i < this->high; i++)
		{
			double otnoshenie = 0;
			for (int k = 0; k < i; k++)
			{
				if (this->value[k][k]!= 0) otnoshenie = this->value[i][k] / this->value[k][k];
				else otnoshenie = 0;
				{
					for (int j = 0; j < this->lengh; j++)
					{
						sub(i, j) -= sub(k, j) * otnoshenie;
						this->value[i][j] -= this->value[k][j] * otnoshenie;
					}
				}
			}
		}
		//�������� ��� �� ��������� �� 0 (��� ��� ��������� �������)
		this->Check_diag_zero(sub);

		//����������� ������ �� ���������
		for (int i = 0; i < this->high; i++)
		{

			for (int j = 0; j < this->lengh; j++)
			{
				if (i != j)
				{
					this->value[i][j] /= this->value[i][i];
					sub(i, j) /= this->value[i][i];
				}
			}
			sub(i, i) /= this->value[i][i];
			this->value[i][i] /= this->value[i][i];
		}

		//����������� � ��������� �������
		for (int i = this->high - 1; i >= 0; i--)
		{
			double otnoshenie = 0;
			for (int k = this->high - 1; k > i; k--)
			{
				otnoshenie = this->value[i][k];
				{
					for (int j = 0; j < this->lengh; j++)
					{
						sub(i, j) -= sub(k, j) * otnoshenie;
						this->value[i][j] -= this->value[k][j] * otnoshenie;
					}
				}
			}
		}
		//���������� ���������
		for (int i = 0; i < this->high; i++)
		{
			for (int j = 0; j < this->lengh; j++)
			{
				this->value[i][j] = sub(i,j);
			}
		}
	}
}

void TMatrix::xoleck()
{
	if (sylver())
	{
		TMatrix sub(this->high, this->lengh);
		//������������� 0
		for (int i = 0; i < this->high; i++)
			for (int j = 0; j < this->high; j++)
				sub(i, j) = 0;
		//������� 0 0
		sub(0, 0) = sqrt(value[0][0]);

		for (int j = 1; j < this->lengh; j++)
			sub(j, 0) = value[0][j] / sub(0, 0);

		for (int i = 1; i < this->lengh; i++)
		{
			double sum = 0;
			for (int p = 0; p < i; p++)
			{
				sum += pow(sub(i, p), 2);
			}
			sub(i, i) = sqrt(value[i][i] - sum);

			if (i != this->lengh - 1)
			{
				double sum2 = 0;
				for (int j = i + 1; j < this->lengh; j++)
				{
					for (int p = 0; p < i; p++)
					{
						sum2 += (sub(i, p)*sub(j, p));
					}
					sub(j, i) = (value[i][j] - sum2) / sub(i, i);
				}
			}
		}
		//��������� ������� �� �����������������
		for (int i = 0; i < this->high; i++)
		{
			for (int j = 0; j < this->lengh; j++)
			{
				this->value[i][j] = sub(i, j);

			}
		}
	}
}

void TMatrix::xoleck_obr()
{
	if (sylver())
	{
		xoleck();
		TMatrix sub = E();

		//��������
		TMatrix Lt(this->high, this->lengh);
		for (int i = 0; i < this->high; i++)
		{
			for (int j = 0; j < this->lengh; j++)
			{
				Lt(i, j) = this->value[j][i];
			}
		}
		TMatrix L(this->high, this->lengh);
		L = Lt;
		//����������
		L = L.transp();

		double sum;
		for (int i = this->lengh - 1; i >= 0; i--)
		{
			for (int j = i; j >= 0; j--)
			{
				sum = 0;
				if (i == j)
				{
					for (int k = i + 1; k < this->lengh; k++)
					{
						double a = sub(k, i);
						sum += L(k, i)*sub(k, i);
					}
					sub(i, j) = ((1 / L(i, j)) - sum) / L(i, j);
				}
				else
				{
					for (int k = j + 1; k < this->lengh; k++)
					{
						double b = sub(k, i);
						sum += L(k, j)*sub(k, i);
					}
					sub(i, j) = -(sum / L(j, j));
					sub(j, i) = sub(i, j);
				}
			}
		}
		//��������� ������� �� ��������
		for (int i = 0; i < this->high; i++)
		{
			for (int j = 0; j < this->lengh; j++)
			{
				this->value[j][i] = sub(i, j);
			}
		}
	}
}

TMatrix TMatrix::E()
{
	TMatrix sub(high, lengh);
	for (int i = 0; i < high; i++)
	{
		sub.value[i] = new double[lengh];
	};
	for (int i = 0; i < high; i++)
	{
		for (int j = 0; j < lengh; j++)
		{
			if (i == j) sub.value[i][j] = 1;
			else sub.value[i][j] = 0;
		}
	}
	return sub;
};

void TMatrix::negative()
{
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			this->value[i][j] = -this->value[i][j];
		}
	}
}
//����������
void TMatrix::resize(int new_high, int new_lengh)
{
	TMatrix sub(new_high, new_lengh);
	int min_high = new_high < this->high ? new_high : this->high;
	int min_lengh = new_lengh < this->lengh ? new_lengh : this->lengh;
	//������ ������������
	for (int i = 0; i < min_high; i++)
	{
		for (int j = 0; j < min_lengh; j++)
		{
			sub.value[i][j] = this->value[i][j];
		}
	}

	if (value)
	{
		for (int i = 0; i < high; i++)
		{
			delete[] value[i];
			value[i] = nullptr;
		}
		delete[] value;
		high = 0;
		lengh = 0;
		value = nullptr;
	}
	this->value = new double*[new_high];
	for (int i = 0; i < new_high; i++)
	{
		this->value[i] = new double[new_lengh];
	}
	high = new_high;
	lengh = new_lengh;
	//�������� ������������
	for (int i = 0; i < min_high; i++)
	{
		for (int j = 0; j < min_lengh; j++)
		{
			this->value[i][j] = sub.value[i][j];
		}
	}
}
// �������� - ������� �����
TMatrix TMatrix::operator - () const
{
	TMatrix sub(high, lengh);
	for (int i = 0; i < this->high; i++)
	{
		for (int j = 0; j < this->lengh; j++)
		{
			sub.value[i][j] = -this->value[i][j];
		}
	}
	return sub;
};
