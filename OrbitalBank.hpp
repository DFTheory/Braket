#include "Elements.hpp"

int Factorial(int x)
{
    double result = 1;
    if (x == 0)
    {
        result = 1;
    }
    else if (x < 0)
    {
        cout << "Error" << endl;
    }
    else
    {
        for (int i = 1; i <= x; i++)
        {
            result *= i;
        }
    }
    return result;
}

double NmlCo(double i, double j, double k, double a) //Calculate the normalization coefficient for basis functions
{
    double PartA = pow(2 * a / M_PI, 0.75);
    double PartB = pow(8 * a, i + j + k) * Factorial(i) * Factorial(j) * Factorial(k);
    double PartC = Factorial(2 * i) * Factorial(2 * j) * Factorial(2 * k);
    double result = PartA * pow(PartB / PartC, 0.5);
    return result;
}

double SGTF(double x, double y, double z, double a)
{
    double N = NmlCo(0, 0, 0, a);
    double result = N * exp(-a * (x * x + y * y + z * z));
    return result;
}

double PxGTF(double x, double y, double z, double a)
{
    double N = NmlCo(1, 0, 0, a);
    double result = N * x * exp(-a * (x * x + y * y + z * z));
    return result;
}

double PyGTF(double x, double y, double z, double a)
{
    double N = NmlCo(0, 1, 0, a);
    double result = N * y * exp(-a * (x * x + y * y + z * z));
    return result;
}

double PzGTF(double x, double y, double z, double a)
{
    double N = NmlCo(0, 0, 1, a);
    double result = N * z * exp(-a * (x * x + y * y + z * z));
    return result;
}

double DxyGTF(double x, double y, double z, double a)
{
    double N = NmlCo(1, 1, 0, a);
    double result = N * x * y * exp(-a * (x * x + y * y + z * z));
    return result;
}

double DxzGTF(double x, double y, double z, double a)
{
    double N = NmlCo(1, 0, 1, a);
    double result = N * x * z * exp(-a * (x * x + y * y + z * z));
    return result;
}

double DyzGTF(double x, double y, double z, double a)
{
    double N = NmlCo(0, 1, 1, a);
    double result = N * y * z * exp(-a * (x * x + y * y + z * z));
    return result;
}

double Dz2GTF(double x, double y, double z, double a)
{
    double N = NmlCo(0, 0, 2, a);
    double result = N * (-0.5 * x * x - 0.5 * y * y + z * z) * exp(-a * (x * x + y * y + z * z));
    return result;
}

double Dx2y2GTF(double x, double y, double z, double a)
{
    double N = NmlCo(2, 0, 0, a);
    double result = N * (sqrt(3) / 2) * (x * x - y * y) * exp(-a * (x * x + y * y + z * z));
    return result;
}