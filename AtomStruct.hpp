#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int PeriodTableHash[118] = {-1};

vector<vector<double>> Orbitals;

struct BSUnit
{
    vector<vector<vector<double>>> S;
    vector<vector<vector<double>>> P;
    vector<vector<vector<double>>> D;
    vector<vector<vector<double>>> F;
};

vector<BSUnit> BSBank;

struct GaussTypeFunction
{
    int xi;
    int yj;
    int zk;
    double ealpha;
    double ecofficient;
};

struct Orbital
{
    vector<GaussTypeFunction> Expression;
    //Add all GTFs inside together u get an expression of the orbital.
};

struct Atom
{
    string AtomSymbol;        //The symbol of an element
    int SerialNumber;         //The sequence number in the input file
    int AtomSerialNumber;     //Atom's sequence number
    vector<double> Cartesian; //In order of xyz.
    vector<double> Sphere;    //In order of r, theta and phi.
    vector<vector<vector<Orbital>>> AtomOrbitals;
    //AtomOrbitals[n-1][m]["l"]. Matches all three quantum numbers, except that u have to pay attention to "l"
};

struct MO
{
    int No;
    double MOenergy;
};