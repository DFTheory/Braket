#include "AtomStruct.hpp"
#include <cmath>
#include <string>
using namespace std;

void CartesianToSphere(vector<Atom> *sample)
{
    for (int i = 0; i < sample->size(); i++)
    {
        double x, y, z;
        x = sample->at(i).Cartesian[0];
        y = sample->at(i).Cartesian[1];
        z = sample->at(i).Cartesian[2];
        double r, theta, phi; // In rad.
        r = sqrt(x * x + y * y + z * z);
        if (r == 0)
        {
            theta = 0;
        }
        else
        {
            theta = acos(z / r);
        }
        if (x == 0)
        {
            if (y > 0)
            {
                phi = M_PI / 2;
            }
            else if (y < 0)
            {
                phi = -M_PI / 2;
            }
            else if (y == 0)
            {
                phi = 0;
            }
        }
        else
        {
            phi = atan(y / x);
        }
        sample->at(i).Sphere.push_back(r);
        sample->at(i).Sphere.push_back(theta);
        sample->at(i).Sphere.push_back(phi);
    }
}
void SphereToCartesian(vector<Atom> *sample)
{
    for (int i = 0; i < sample->size(); i++)
    {
        if (sample->at(i).Sphere.size() == 0)
        {
            cout << "Error when transforming from Sphere to Cartesian at atom " << sample->at(i).SerialNumber << endl;
        }
        else
        {
            double r, theta, phi;
            r = sample->at(i).Sphere[0];
            theta = sample->at(i).Sphere[1];
            phi = sample->at(i).Sphere[2];
            sample->at(i).Cartesian[0] = r * sin(theta) * cos(phi);
            sample->at(i).Cartesian[1] = r * sin(theta) * sin(phi);
            sample->at(i).Cartesian[2] = r * cos(theta);
        }
    }
}