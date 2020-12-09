#pragma once
#include "VectMatr.h"
#include <vector>
#include <iostream>
#include <math.h>

//class TQuaternion;

class TVector
{
	double* value;
	int size;
public:
	//����������� �� ���������
	TVector(){ value = nullptr; size = 0; };
	//����������� �������������� �������
	TVector(int);
	//����������
	~TVector();
	// ����������� �����������
	TVector(const TVector& rvalue);
	//������������
	void operator = (const TVector&);
	//������� ���������� ��������� � ������ ����������� �������
	void add(double []);
	//������� ���������� ��������� � ������ ����������� �������
	void add(std::vector<double> vec);
	//������� ��������� �������� ����������� ��������
	double& operator [] (int index) const { return value[index]; };
	//������� ��������� �������� ������� �������
	int length() const { return size; };
	//������� ��������� ������
	double mod() const;
	//������� ��������� �������������� �������
	TVector norm() const;
	//������� ��������� �����
	TVector negative();
	//�������� ��������
	TVector operator + (const TVector&) const;
	//��������� ��������
	TVector operator - (const TVector&) const;
	//��������� ������������
	double operator * (const TVector&) const;
	//��������� �� �����
	TVector operator * (const double&) const;
	//��������� ������������
	TVector operator ^ (const TVector&) const;
	//����������
	void resize(int new_size);
	// ������������� ������� - �������� ��������� ����� �� ������
	friend TVector operator * (double lvalue, const TVector& rvalue){ { TVector sub(rvalue.size); sub = rvalue; return sub * lvalue; }; };
	TVector operator - () const;
	// ������� ������� ������ �������� ��� �� �������� ���� ��� ������ ������� �������
	TVector Rodrig(double phi, const TVector& axis);
};
class TMatrix
{
	double **value;
	int high, lengh;
public:
	//����������� �� ���������
	TMatrix(){ high = 0; lengh = 0; };
	//����������� �������������� �������
	TMatrix(int, int);
	//����������
	~TMatrix();
	// ����������� �����������
	TMatrix(const TMatrix& rvalue);
	//������� ������ � ������
	int getHigh(){ return high; };
	int getLengh(){ return lengh; };
	//������� ��������� �������� ����������� ��������
	double& operator () (int index1, int index2) const { return this->value[index1][index2]; };
	//������������
	void operator = (const TMatrix&);
	TMatrix operator + (const TMatrix&);
	TMatrix operator - (const TMatrix&);
	TMatrix operator * (const TMatrix&);
	TVector operator * (const TVector&);
	TMatrix operator * (const double&);                                  

	int search(double, bool, unsigned int&, unsigned int&, unsigned int, unsigned int);
	//����� ��������
	void swapcolumns(unsigned int, unsigned int);
	//����� �����
	void swaprows(unsigned int, unsigned int);
	//�����������
	double det();
	//�������� ����������
	bool sylver();
	//������� ���������� ��������� � ������� ����������� �������
	void add(double*[]);
	//������� ���������� ��������� � ������� ����������� ��������
	void add(std::vector<double> vec, int length);
	TMatrix transp();
	void Check_diag_zero(TMatrix& sub);
	void gaus();
	void xoleck();
	void xoleck_obr();
	TMatrix E();
	//������� ��������� �����
	void negative();
	//����������
	void resize(int new_high, int new_lenght);
	// ������������� ������� - �������� ��������� ����� �� �������
	friend TMatrix operator * (double lvalue, TMatrix& rvalue){ TMatrix sub(rvalue.high, rvalue.lengh); sub = rvalue; return sub * lvalue; };
	// �������� - ������� �����
	TMatrix operator - () const;
};