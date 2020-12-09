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
	//конструктор по умолчанию
	TVector(){ value = nullptr; size = 0; };
	//конструктор фиксированного размера
	TVector(int);
	//деструктор
	~TVector();
	// Конструктор копирования
	TVector(const TVector& rvalue);
	//присваивание
	void operator = (const TVector&);
	//функция добавления элементов в вектор посредством массива
	void add(double []);
	//функция добавления элементов в вектор посредством вектора
	void add(std::vector<double> vec);
	//функция получения значения конкретного элемента
	double& operator [] (int index) const { return value[index]; };
	//функция получения текущего размера вектора
	int length() const { return size; };
	//функция получения модуля
	double mod() const;
	//функция получения нормированного вектора
	TVector norm() const;
	//функция обращения знака
	TVector negative();
	//сложение векторов
	TVector operator + (const TVector&) const;
	//вычитание векторов
	TVector operator - (const TVector&) const;
	//скалярное произведение
	double operator * (const TVector&) const;
	//умножение на число
	TVector operator * (const double&) const;
	//векторное произведение
	TVector operator ^ (const TVector&) const;
	//переразмер
	void resize(int new_size);
	// Дружественная функция - оператор умножения числа на вектор
	friend TVector operator * (double lvalue, const TVector& rvalue){ { TVector sub(rvalue.size); sub = rvalue; return sub * lvalue; }; };
	TVector operator - () const;
	// Поворот вектора вокруг заданной оси на заданный угол при помощи формулы Родрига
	TVector Rodrig(double phi, const TVector& axis);
};
class TMatrix
{
	double **value;
	int high, lengh;
public:
	//конструктор по умолчанию
	TMatrix(){ high = 0; lengh = 0; };
	//конструктор фиксированного размера
	TMatrix(int, int);
	//деструктор
	~TMatrix();
	// Конструктор копирования
	TMatrix(const TMatrix& rvalue);
	//возврат высоты и длинны
	int getHigh(){ return high; };
	int getLengh(){ return lengh; };
	//функция получения значения конкретного элемента
	double& operator () (int index1, int index2) const { return this->value[index1][index2]; };
	//присваивание
	void operator = (const TMatrix&);
	TMatrix operator + (const TMatrix&);
	TMatrix operator - (const TMatrix&);
	TMatrix operator * (const TMatrix&);
	TVector operator * (const TVector&);
	TMatrix operator * (const double&);                                  

	int search(double, bool, unsigned int&, unsigned int&, unsigned int, unsigned int);
	//смена столбцов
	void swapcolumns(unsigned int, unsigned int);
	//смена строк
	void swaprows(unsigned int, unsigned int);
	//детерминант
	double det();
	//критерий сильвестра
	bool sylver();
	//функция добавления элементов в матрицу посредством массива
	void add(double*[]);
	//функция добавления элементов в матрицу посредством вектором
	void add(std::vector<double> vec, int length);
	TMatrix transp();
	void Check_diag_zero(TMatrix& sub);
	void gaus();
	void xoleck();
	void xoleck_obr();
	TMatrix E();
	//функция обращения знака
	void negative();
	//переразмер
	void resize(int new_high, int new_lenght);
	// Дружественная функция - оператор умножения числа на матрицу
	friend TMatrix operator * (double lvalue, TMatrix& rvalue){ TMatrix sub(rvalue.high, rvalue.lengh); sub = rvalue; return sub * lvalue; };
	// Оператор - унарный минус
	TMatrix operator - () const;
};