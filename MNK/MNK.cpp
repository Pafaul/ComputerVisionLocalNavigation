#include <iostream>
#include "VectMatr.h"
#include <conio.h>
#include <string>
#include <fstream>

float val_rand = 0;

void testMNK_final(std::vector<float> x, std::vector<float> y, std::vector<float> z, std::vector<float> X, std::vector<float> Y, std::vector<float> Z)
{
	TMatrix R(3 * x.size(), 12);

	std::vector<double> array_r;
	for (int i = 0; i < x.size(); i++)
	{
		array_r.push_back(x[i]);//r11
		array_r.push_back(y[i]);//r12
		array_r.push_back(z[i]);//r13
		array_r.push_back(0);//r21
		array_r.push_back(0);//r22
		array_r.push_back(0);//r23
		array_r.push_back(0);//r31
		array_r.push_back(0);//r32
		array_r.push_back(0);//r33
		array_r.push_back(1);//hx
		array_r.push_back(0);//hy
		array_r.push_back(0);//hz

		array_r.push_back(0);//r11
		array_r.push_back(0);//r12
		array_r.push_back(0);//r13
		array_r.push_back(x[i]);//r21
		array_r.push_back(y[i]);//r22
		array_r.push_back(z[i]);//r23
		array_r.push_back(0);//r31
		array_r.push_back(0);//r32
		array_r.push_back(0);//r33
		array_r.push_back(0);//hx
		array_r.push_back(1);//hy
		array_r.push_back(0);//hz

		array_r.push_back(0);//r11
		array_r.push_back(0);//r12
		array_r.push_back(0);//r13
		array_r.push_back(0);//r21
		array_r.push_back(0);//r22
		array_r.push_back(0);//r23
		array_r.push_back(x[i]);//r31
		array_r.push_back(y[i]);//r32
		array_r.push_back(z[i]);//r33
		array_r.push_back(0);//hx
		array_r.push_back(0);//hy
		array_r.push_back(1);//hz
	}
	R.add(array_r, 12);

	TVector vX(3 * X.size());
	std::vector<double> vXX;
	
	for (int i = 0; i < X.size(); i++)
	{
		vXX.push_back(X[i]);
		vXX.push_back(Y[i]);
		vXX.push_back(Z[i]);
	}

	vX.add(vXX);

	TVector Res(12);
	//R.transp();
	TMatrix Rtrans(2, 2);
	Rtrans = R.transp();
	//(R.transp() * R);
	TMatrix temp(2, 2);
	temp = (Rtrans*R);
	//draw(temp);
	//(R.transp() * R).gaus();
	temp.gaus();
	//(R.transp() * R).gaus() * R.transp();
	temp = temp * Rtrans;
	//(R.transp() * R).gaus() * R.transp() * X;
	Res = temp * vX;

	float alpha = asin(Res[2]);
	float gamma = 0;
	float omega = 0;
	if (abs(cos(alpha)) >= 1e-8)
	{
		if (fabs(Res[0] / cos(alpha) - 1) < 1e-7) gamma = 0;
		else gamma = acos(Res[0] / cos(alpha));
		if (fabs(Res[8] / cos(alpha) - 1) < 1e-7) omega = 0;
		else omega = acos(Res[8] / cos(alpha));
	}
	else
	{
		gamma = 0;
		omega = 0;
		std::cout << "alpha = 0, gamma and omega havenot answer";
	}

	std::cout << "alpha = " << std::to_string(alpha * 180 / 3.1416) << " gamma = " << std::to_string(gamma * 180 / 3.1416) << " omega = " << std::to_string(omega * 180 / 3.1416) << std::endl;
	std::cout << "hx = " << Res[9] << " hy = " << Res[10] << " hz = " << Res[11] << std::endl;

	std::ofstream out("angle_and_translation.txt", std::ios::app);
	//out.precision(2);
	out << alpha * 180 / 3.1416f << " " << gamma * 180 / 3.1416f << " " << omega * 180 / 3.1416f << " " 
		<< Res[9] << " " << Res[10] << " " << Res[11] << std::endl;
	out.close();
}

void create_point(std::vector<float>& vx, std::vector<float>& vy, std::vector<float>& vz, std::vector<float>& vX, std::vector<float>& vY, std::vector<float>& vZ)
{
	float a = 0;
	float g = 0;
	float w = -8;
	a = a / 180 * 3.1416;
	g = g / 180 * 3.1416;
	w = w / 180 * 3.1416;

	float hx = 20;
	float hy = 10;
	float hz = 0;

	for (int i = 0; i < 20; i++)
	{
		int randx = std::rand() % 200;
		int randy = std::rand() % 200;
		int randz = std::rand() % 200;

		float x = randx;
		float y = randy;
		float z = randz;

		float ex = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand
		float ey = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand
		float ez = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand

		vx.push_back(x + ex);
		vy.push_back(y + ey);
		vz.push_back(z + ez);

		TVector X(3);
		TMatrix A(3, 3);
		double m_r1[] = { cos(a)*cos(g), -cos(a)*sin(g), sin(a) };
		double m_r2[] = { sin(w)*sin(a)*cos(g)+cos(w)*sin(g), -sin(w)*sin(a)*sin(g)+cos(w)*cos(g), -sin(w)*cos(a) };
		double m_r3[] = { -cos(w)*sin(a)*cos(g), cos(w)*sin(a)*sin(g)+sin(w)*cos(g), cos(w)*cos(a) };
		double* m_r[] = { m_r1, m_r2, m_r3 };
		A.add(m_r);

		TVector Tvx(3);
		double vXX[] = { x, y, z };
		Tvx.add(vXX);

		X = A * Tvx;

		float ex2 = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand
		float ey2 = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand
		float ez2 = (std::rand() % (200) - 100) / 100.f * val_rand;//-val_rand:val_rand

		vX.push_back(X[0] + hx + ex2);
		vY.push_back(X[1] + hy + ey2);
		vZ.push_back(X[2] + hz + ez2);
	}

	std::ofstream out("input.txt");
	for (int i = 0; i < vx.size(); i++)
	{
		out << vx[i] << " " << vy[i] << " " << vz[i] << " " << vX[i] << " " << vY[i] << " " << vZ[i] << std::endl;
	}
	out.close();
}

void set_point(std::vector<float>& vx, std::vector<float>& vy, std::vector<float>& vz, std::vector<float>& vX, std::vector<float>& vY, std::vector<float>& vZ, std::string in)
{
	std::ifstream input(in);
	float x, y, z, X, Y, Z;
	while (input >> x >> y >> z >> X >> Y >> Z)
	{
		vx.push_back(x);
		vy.push_back(y);
		vz.push_back(z);
		vX.push_back(X);
		vY.push_back(Y);
		vZ.push_back(Z);
	}

}

void main(int argc, char* argv[])
{

	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> z;
	std::vector<float> X;
	std::vector<float> Y;
	std::vector<float> Z;

	std::string input = "input.txt";
	if (argc > 1) input = argv[1];

	std::ofstream out("angle_and_translation.txt");
	out.close();

	for (int i = 0; i < 50; i++)
	{
		val_rand = 0;
		x.clear();
		y.clear();
		z.clear();
		X.clear();
		Y.clear();
		Z.clear();
		create_point(x, y, z, X, Y, Z);
		//set_point(x, y, z, X, Y, Z, input);
		testMNK_final(x, y, z, X, Y, Z);
	}
	for (int i = 0; i < 50; i++)
	{
		val_rand = sqrt(3);
		x.clear();
		y.clear();
		z.clear();
		X.clear();
		Y.clear();
		Z.clear();

		create_point(x, y, z, X, Y, Z);
		//set_point(x, y, z, X, Y, Z, input);
		testMNK_final(x, y, z, X, Y, Z);
	}
	for (int i = 0; i < 50; i++)
	{
		val_rand = sqrt(6);
		x.clear();
		y.clear();
		z.clear();
		X.clear();
		Y.clear();
		Z.clear(); 
		
		create_point(x, y, z, X, Y, Z);
		//set_point(x, y, z, X, Y, Z, input);
		testMNK_final(x, y, z, X, Y, Z);
	}
	for (int i = 0; i < 50; i++)
	{
		val_rand = 3;
		x.clear();
		y.clear();
		z.clear();
		X.clear();
		Y.clear();
		Z.clear(); 
		
		create_point(x, y, z, X, Y, Z);
		//set_point(x, y, z, X, Y, Z, input);
		testMNK_final(x, y, z, X, Y, Z);
	}
	//cv::waitKey(0);
}