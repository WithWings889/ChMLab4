#include <iostream>
#include <vector>
#define PI 3.14159265358979323846
#include <locale>
#include "windows.h"

using namespace std;

double F(double x)
{
	return 2 / x;
}
double exaltationInDegree(double a, int x)
{
	double ans = 1.0;
	for (int i = 1; i <= x; ++i)
	{
		ans *= a;
	}
	return ans;
}

int Factorial(int n)
{
	int ans = 1;
	for (int i = 1; i <= n; ++i)
	{
		ans *= i;
	}
	return ans;
}

double M11(double a, double b)
{
	if (a > b)
		return abs(-2 * (double)Factorial(11) / exaltationInDegree(b, 12));
	else
		return abs(-2 * (double)Factorial(11) / exaltationInDegree(a, 12));
}

void fildWithDataOnVectorChebysh(double a, double b, int k, vector<double>& x, vector<double>& y)
{
	for (int i = 0; i < k; ++i)
	{
		x[i] =(a + b)/2 + (b - a)/2 * cos((2*i + 1)*PI/(2*k));
		y[i] = F(x[i]);
	}
}

void fildWithDataOnVectorEquidistant(double a, double b, int k, vector<double>& x, vector<double>& y)
{
	double h = (b - a) / k;
	for (int i = 0; i < k; ++i)
	{
		x[i] = a + h * i;
		y[i] = F(x[i]);
	}
}

double Lagrange(vector<double> x, vector<double> y, int n, double _x)
{
	double result = 0.0;

	for (short i = 0; i < n; i++)
	{
		double P = 1.0;

		for (short j = 0; j < n; j++)
			if (j != i)
				P *= (_x - x[j]) / (x[i] - x[j]);

		result += P * y[i];
	}

	return result;
}


double accuracyChebysh(double a, double b, int k)
{
	return M11(a, b) /Factorial(k + 1) * exaltationInDegree(b - a, k + 1) / exaltationInDegree(2, 2 * k + 1);
}

double accuracyEquidistant(double a, double b, int k, double _x, vector<double> x)
{
	double w = 1;
	for (int i = 0; i < k; ++i)
	{
		w *= abs(_x - x[i]);
	}
	return M11(a, b) /Factorial(k + 1) * w;
}

void AnsLagrangeChebysh(double a, double b, int k, double _x)
{
	vector<double> x(k), y(k);
	fildWithDataOnVectorChebysh(a, b, k, x, y);
	cout << "Поліном Лагранджа за нулями Чебишова "  << endl;
	cout << "Значення в точці " << Lagrange(x, y, k, _x) << endl;
	cout << "Точність " << accuracyChebysh(a, b, k) << endl;
	cout << endl;
}

void AnsLagrangeEquidistant(double a, double b, int k, double _x)
{
	vector<double> p(k), q(k);
	fildWithDataOnVectorEquidistant(a, b, k, p, q);
	cout << "Поліном Лагранджа за рівновіддаленили вузлами " << endl;
	cout << "Значення в точці " << Lagrange(p, q, k, _x) << endl;
	cout << "Точність " << accuracyEquidistant(a, b, k, _x, p) << endl;
	cout << endl;
}

void printFunc(double _x)
{
	cout << "Оберемо функцію  на проміжку в точці 2/х"  << endl;
	cout << "Значення функції в точці " << F(_x) << endl;
	cout << endl;
}

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	double a = 1, b = 3, _x = 2.0876;
	int k = 10;
	printFunc(_x);
	AnsLagrangeChebysh(a, b, k, _x);
	AnsLagrangeEquidistant(a, b, k, _x);
}