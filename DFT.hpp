#include "Functionals.hpp"
using namespace std;
using namespace Eigen;

double FCTETensor[100][100][100][100];

//Step 1: Calculate the repulsive energy of nucleus
//Step 2: Calculate the overlap integral matrix
//Step 3: Calculate the kinetic energy integral matrix
//Step 4: Calculate the nuclear attraction energy matrix
//Step 5: Calculate the 4 center 2 electron integrals tensor
//Step 6: Prepare for the self-consistent field procedure
//Step 7: Run the self-consistent calculation until converge
//Step 8: Output results
void DFT(double Energy, vector<Atom> inp, int Electrons, string FileName, int FuntionalType)
{
    fstream DFTout(FileName + ".out", ios::out);
    double RepulsiveEnergy = NuclearRepulsiveEnergy(inp);
    Matrix<double, Dynamic, Dynamic> OverlapIntegralMatrix;
    int NumberOfOrbitals = CountOrbitals(inp);
    OverlapIntegralMatrix.resize(NumberOfOrbitals, NumberOfOrbitals);
    QuantumNumber OrbitalSequence[100];
    //Create orbital sequence
    int now = 0;
    for (int i = 0; i < inp.size(); i++)
    {
        for (int j = 0; j < inp[i].AtomOrbitals.size(); j++)
        {
            for (int k = 0; k < inp[i].AtomOrbitals[j].size(); k++)
            {
                for (int l = 0; l < inp[i].AtomOrbitals[j][k].size(); l++)
                {
                    OrbitalSequence[now].SerialNumber = inp[i].SerialNumber;
                    OrbitalSequence[now].ProtonNumber = inp[i].AtomSerialNumber;
                    OrbitalSequence[now].n = j + 1;
                    OrbitalSequence[now].l = k;
                    OrbitalSequence[now].m = l - k;
                    now++;
                }
            }
        }
    }
    //Calculate overlap matrix
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = i; j < NumberOfOrbitals; j++)
        {
            int A = OrbitalSequence[i].SerialNumber - 1;
            int B = OrbitalSequence[j].SerialNumber - 1;
            int nA = OrbitalSequence[i].n - 1;
            int nB = OrbitalSequence[j].n - 1;
            int lA = OrbitalSequence[i].l;
            int lB = OrbitalSequence[j].l;
            int mA = OrbitalSequence[i].m + lA;
            int mB = OrbitalSequence[j].m + lB;
            OverlapIntegralMatrix(i, j) = OrbitalSIntegral(inp[A], inp[B], inp[A].AtomOrbitals[nA][lA][mA], inp[B].AtomOrbitals[nB][lB][mB]);
            if (OverlapIntegralMatrix(i, j) < pow(10, -9) && OverlapIntegralMatrix(i, j) > -pow(10, -9))
            {
                OverlapIntegralMatrix(i, j) = 0;
            }
            if (i != j)
            {
                OverlapIntegralMatrix(j, i) = OverlapIntegralMatrix(i, j);
            }
        }
    }
    //Calculate kinetic energy matrix
    Matrix<double, Dynamic, Dynamic> KineticEnergyMatrix;
    KineticEnergyMatrix.resize(NumberOfOrbitals, NumberOfOrbitals);
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = i; j < NumberOfOrbitals; j++)
        {
            int A = OrbitalSequence[i].SerialNumber - 1;
            int B = OrbitalSequence[j].SerialNumber - 1;
            int nA = OrbitalSequence[i].n - 1;
            int nB = OrbitalSequence[j].n - 1;
            int lA = OrbitalSequence[i].l;
            int lB = OrbitalSequence[j].l;
            int mA = OrbitalSequence[i].m + lA;
            int mB = OrbitalSequence[j].m + lB;
            KineticEnergyMatrix(i, j) = OrbitalKineticEnergyIntegral(inp[A], inp[B], inp[A].AtomOrbitals[nA][lA][mA], inp[B].AtomOrbitals[nB][lB][mB]);
            if (KineticEnergyMatrix(i, j) < pow(10, -9) && KineticEnergyMatrix(i, j) > -pow(10, -9))
            {
                KineticEnergyMatrix(i, j) = 0;
            }
            if (i != j)
            {
                KineticEnergyMatrix(j, i) = KineticEnergyMatrix(i, j);
            }
        }
    }
    //Calculate nuclear attraction energy matrix
    Matrix<double, Dynamic, Dynamic> NuclearAttractionMatrix;
    NuclearAttractionMatrix.resize(NumberOfOrbitals, NumberOfOrbitals);
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = i; j < NumberOfOrbitals; j++)
        {
            int A = OrbitalSequence[i].SerialNumber - 1;
            int B = OrbitalSequence[j].SerialNumber - 1;
            int nA = OrbitalSequence[i].n - 1;
            int nB = OrbitalSequence[j].n - 1;
            int lA = OrbitalSequence[i].l;
            int lB = OrbitalSequence[j].l;
            int mA = OrbitalSequence[i].m + lA;
            int mB = OrbitalSequence[j].m + lB;
            NuclearAttractionMatrix(i, j) = OrbitalNuclearAttractionIntegral(inp[A], inp[B], inp[A].AtomOrbitals[nA][lA][mA], inp[B].AtomOrbitals[nB][lB][mB], inp);
            if (NuclearAttractionMatrix(i, j) < pow(10, -9) && NuclearAttractionMatrix(i, j) > -pow(10, -9))
            {
                NuclearAttractionMatrix(i, j) = 0;
            }
            if (i != j)
            {
                NuclearAttractionMatrix(j, i) = NuclearAttractionMatrix(i, j);
            }
        }
    }
    //Calculate the 4 center 2 electron integrals
    //Calculate the 4 center 2 electron integrals
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            for (int p = 0; p < 100; p++)
            {
                for (int q = 0; q < 100; q++)
                {
                    FCTETensor[i][j][p][q] = 100;
                }
            }
        }
    }
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = 0; j < NumberOfOrbitals; j++)
        {
            for (int p = 0; p < NumberOfOrbitals; p++)
            {
                for (int q = 0; q < NumberOfOrbitals; q++)
                {
                    if (FCTETensor[i][j][p][q] == 100)
                    {
                        int A, B, C, D;
                        int nA, nB, nC, nD;
                        int lA, lB, lC, lD;
                        int mA, mB, mC, mD;
                        A = OrbitalSequence[i].SerialNumber - 1;
                        B = OrbitalSequence[j].SerialNumber - 1;
                        C = OrbitalSequence[p].SerialNumber - 1;
                        D = OrbitalSequence[q].SerialNumber - 1;
                        nA = OrbitalSequence[i].n - 1;
                        nB = OrbitalSequence[j].n - 1;
                        nC = OrbitalSequence[p].n - 1;
                        nD = OrbitalSequence[q].n - 1;
                        lA = OrbitalSequence[i].l;
                        lB = OrbitalSequence[j].l;
                        lC = OrbitalSequence[p].l;
                        lD = OrbitalSequence[q].l;
                        mA = OrbitalSequence[i].m + lA;
                        mB = OrbitalSequence[j].m + lB;
                        mC = OrbitalSequence[p].m + lC;
                        mD = OrbitalSequence[q].m + lD;
                        double MatrixElement = OrbitalFCTEIntegral(inp[A], inp[B], inp[C], inp[D], inp[A].AtomOrbitals[nA][lA][mA], inp[B].AtomOrbitals[nB][lB][mB], inp[C].AtomOrbitals[nC][lC][mC], inp[D].AtomOrbitals[nD][lD][mD]);
                        if (MatrixElement < pow(10, -9) && MatrixElement > -pow(10, -9))
                        {
                            MatrixElement = 0;
                        }
                        FCTETensor[i][j][p][q] = MatrixElement;
                        FCTETensor[j][i][p][q] = MatrixElement;
                        FCTETensor[i][j][q][p] = MatrixElement;
                        FCTETensor[j][i][q][p] = MatrixElement;
                        FCTETensor[p][q][i][j] = MatrixElement;
                        FCTETensor[p][q][j][i] = MatrixElement;
                        FCTETensor[q][p][i][j] = MatrixElement;
                        FCTETensor[q][p][j][i] = MatrixElement;
                    }
                }
            }
        }
    }
    //SCF Preparation
    //Diagonalize overlap matrix
}