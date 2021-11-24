#include <omp.h>
#include "Run.hpp"
using namespace std;

int main(int argc, char **argv)
{
    //Parameters of the system.
    vector<Atom> Molecule;
    int SpinMultiplicity, Charge, NumberOfAtoms;
    int ThreadNumber;
    int JobType = 0, BasisSet = -1, Method = 0;
    int NumberOfElectrons = 0;
    double TotalEnergy = 0;
    string InputFileName = argv[1];
    string TheFileName = InputFileName.substr(0, InputFileName.size() - 3);
    fstream fin(argv[1], ios::in);
    string jt, mtd, bs;
    fin >> jt >> mtd >> bs;
    fin >> Charge >> SpinMultiplicity;
    fin >> NumberOfAtoms;
    for (int i = 0; i < NumberOfAtoms; i++)
    {
        Atom temp;
        string element;
        double x, y, z;
        fin >> element >> temp.SerialNumber >> x >> y >> z;
        temp.Cartesian.push_back(x * 1.8897258);
        temp.Cartesian.push_back(y * 1.8897258);
        temp.Cartesian.push_back(z * 1.8897258);
        temp.AtomSerialNumber = AtomSelector(element);
        temp.AtomSymbol = element;
        Molecule.push_back(temp);
        NumberOfElectrons += temp.AtomSerialNumber;
    }
    NumberOfElectrons -= Charge;
    JobType = JobControl(jt);
    Method = MethodControl(mtd);
    BasisSet = BasisSetsControl(bs);
    if (JobType == 0 || Method == 0 || BasisSet == -1 || NumberOfAtoms == 0)
    {
        cout << "Please modify the input file and retry";
    }
    else
    {
        BuildBSBank(BSBank, BasisSet);
        for (int i = 0; i < Molecule.size(); i++)
        {
            OrbitalMaker(BasisSet, Molecule[i].AtomSerialNumber, Molecule[i]);
        }
        //HartreeFock(TotalEnergy, Molecule, NumberOfElectrons, TheFileName);
        RunMission(JobType, Method, BasisSet, TotalEnergy, Molecule, NumberOfElectrons, TheFileName);
    }
    return 0;
}