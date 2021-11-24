#include "MathBank.hpp"
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <algorithm>
#include <cmath>
using namespace std;
using namespace Eigen;

double gabcd[100][100][100][100];

struct QuantumNumber
{
    int ProtonNumber;
    int SerialNumber;
    int n, l, m; //m is slightly different for the case of d orbitals.
};

int CountOrbitals(vector<Atom> inp)
{
    int result = 0;
    for (int a = 0; a < inp.size(); a++)
    {
        for (int b = 0; b < inp[a].AtomOrbitals.size(); b++)
        {
            for (int c = 0; c < inp[a].AtomOrbitals[b].size(); c++)
            {
                result += inp[a].AtomOrbitals[b][c].size();
            }
        }
    }
    return result;
}

double NuclearRepulsiveEnergy(vector<Atom> inp)
{
    double result = 0;
    for (int i = 0; i < inp.size(); i++)
    {
        for (int j = i + 1; j < inp.size(); j++)
        {
            double r = 0;
            for (int k = 0; k < 3; k++)
            {
                r += pow(inp[i].Cartesian[k] - inp[j].Cartesian[k], 2);
            }
            r = sqrt(r);
            result += inp[i].AtomSerialNumber * inp[j].AtomSerialNumber / r;
        }
    }
    return result;
}

void NormalizationC(Matrix<double, Dynamic, Dynamic> &Coe, Matrix<double, Dynamic, Dynamic> SMatrix)
{
    for (int i = 0; i < Coe.cols(); i++)
    {
        double SumOfIt = 0;
        for (int j = 0; j < Coe.rows(); j++)
        {
            for (int k = 0; k < Coe.rows(); k++)
            {
                SumOfIt += Coe(k, i) * Coe(j, i) * SMatrix(k, j);
            }
        }
        double NCoe = sqrt(1 / SumOfIt);
        for (int j = 0; j < Coe.rows(); j++)
        {
            Coe(j, i) *= NCoe;
        }
    }
}

void PrintOrbital(int SN, int n, int l, int m, string &TheOrbital, vector<Atom> InputMolecule)
{
    TheOrbital = to_string(SN);
    TheOrbital += InputMolecule[SN - 1].AtomSymbol;
    TheOrbital += to_string(n);
    string OrbitalType;
    if (l == 0)
    {
        OrbitalType = "s";
    }
    else if (l == 1)
    {
        OrbitalType = "p";
        if (m == -1)
        {
            OrbitalType += "x";
        }
        else if (m == 0)
        {
            OrbitalType += "y";
        }
        else if (m == 1)
        {
            OrbitalType += "z";
        }
    }
    else if (l == 2)
    {
        OrbitalType = "d";
        if (m == -2)
        {
            OrbitalType += "dz2";
        }
        else if (m == -1)
        {
            OrbitalType += "dxz";
        }
        else if (m == 0)
        {
            OrbitalType += "dyz";
        }
        else if (m == 1)
        {
            OrbitalType += "dx2-y2";
        }
        else if (m == 2)
        {
            OrbitalType += "dxy";
        }
    }
    TheOrbital += OrbitalType;
}

bool CompareEnergy(pair<int, double> a, pair<int, double> b)
{
    return a.second < b.second;
}

void ArrangeMOs(vector<pair<int, double>> arranged, Matrix<double, Dynamic, Dynamic> in, Matrix<double, Dynamic, Dynamic> &out, int NOO)
{
    for (int i = 0; i < NOO; i++)
    {
        for (int j = 0; j < NOO; j++)
        {
            out(j, i) = in(j, arranged[i].first);
        }
    }
}

void PrintMatrixForMatlab(Matrix<double, Dynamic, Dynamic> A, string name)
{
    cout << name << "=[";
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            cout << A(i, j);
            if (j < A.cols() - 1)
            {
                cout << " ";
            }
            else
            {
                if (i < A.rows() - 1)
                {
                    cout << ";";
                }
            }
        }
        if (i < A.rows() - 1)
        {
            cout << endl;
        }
        else
        {
            cout << "]" << endl;
        }
    }
}

//Step 1: Calculate the repulsive energy of nucleus
//Step 2: Calculate the overlap integral matrix
//Step 3: Calculate the kinetic energy integral matrix and nuclear attraction energy matrix
//Step 4: Calculate the 4 center 2 electron integrals
//Step 5: Preparation for Self-consistent field calculation, including the S^(-0.5) matrix, P matrix
//Step 6: SCF calculation until energy converges

void HartreeFock(double Energy, vector<Atom> inp, int Electrons, string FileName)
{
    fstream HFout(FileName + ".out", ios::out);
    double RepulsiveEnergy = NuclearRepulsiveEnergy(inp);
    Matrix<double, Dynamic, Dynamic> OverlapIntegralMatrix;
    int NumberOfOrbitals = CountOrbitals(inp);
    OverlapIntegralMatrix.resize(NumberOfOrbitals, NumberOfOrbitals);
    QuantumNumber OrbitalSequence[100];
    //Create orbital sequence
    int nowx = 0, nowy = 0, now = 0;
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
    //Overlap matrix
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
    //Kinetic energy matrix
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
    //Nuclear attraction energy matrix
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
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            for (int p = 0; p < 100; p++)
            {
                for (int q = 0; q < 100; q++)
                {
                    gabcd[i][j][p][q] = 100;
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
                    if (gabcd[i][j][p][q] == 100)
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
                        gabcd[i][j][p][q] = MatrixElement;
                        gabcd[j][i][p][q] = MatrixElement;
                        gabcd[i][j][q][p] = MatrixElement;
                        gabcd[j][i][q][p] = MatrixElement;
                        gabcd[p][q][i][j] = MatrixElement;
                        gabcd[p][q][j][i] = MatrixElement;
                        gabcd[q][p][i][j] = MatrixElement;
                        gabcd[q][p][j][i] = MatrixElement;
                    }
                }
            }
        }
    }
    //SCF preparation
    //Diagonalize S matrix
    Matrix<double, Dynamic, Dynamic> Sminus, Spostive, SEValminus, SEValpostive, SEigenVectors, SEigenValues, SEVectorsTranspose; //The S^(-0.5) matrix and S^(0.5) matrix.
    SEigenValues.resize(NumberOfOrbitals, NumberOfOrbitals);
    SEValminus.resize(NumberOfOrbitals, NumberOfOrbitals);
    SEValpostive.resize(NumberOfOrbitals, NumberOfOrbitals);
    EigenSolver<MatrixXd> DiagonalizeS(OverlapIntegralMatrix);
    SEigenVectors = DiagonalizeS.eigenvectors().real();
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = 0; j < NumberOfOrbitals; j++)
        {
            if (i != j)
            {
                SEigenValues(i, j) = 0;
                SEValpostive(i, j) = 0;
                SEValminus(i, j) = 0;
            }
            else
            {
                SEigenValues(i, i) = DiagonalizeS.eigenvalues()[i].real();
                SEValpostive(i, i) = sqrt(SEigenValues(i, i));
                SEValminus(i, i) = 1 / SEValpostive(i, i);
            }
        }
    }
    SEVectorsTranspose = SEigenVectors.transpose();
    Sminus = SEigenVectors * SEValminus * SEVectorsTranspose;
    Spostive = SEigenVectors * SEValpostive * SEVectorsTranspose;
    //Initialize the C matrix
    Matrix<double, Dynamic, Dynamic> CMatrix;
    Matrix<double, Dynamic, Dynamic> CMatrixTTemp;
    Matrix<double, Dynamic, Dynamic> HCore;
    HCore = KineticEnergyMatrix + NuclearAttractionMatrix;
    Matrix<double, Dynamic, Dynamic> Fock;
    Matrix<double, Dynamic, Dynamic> FockT, CMatrixT;
    FockT.resize(NumberOfOrbitals, NumberOfOrbitals);
    CMatrixT.resize(NumberOfOrbitals, NumberOfOrbitals);
    Fock = HCore;
    EigenSolver<MatrixXd> DiagonalizeFock(Sminus * Fock * Sminus);
    CMatrix.resize(NumberOfOrbitals, NumberOfOrbitals);
    CMatrixTTemp.resize(NumberOfOrbitals, NumberOfOrbitals);
    vector<pair<int, double>> MOEnergy;
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        pair<int, double> TEMP(i, DiagonalizeFock.eigenvalues()[i].real());
        MOEnergy.push_back(TEMP);
    }
    CMatrixTTemp = DiagonalizeFock.eigenvectors().real();
    stable_sort(MOEnergy.begin(), MOEnergy.end(), CompareEnergy);
    ArrangeMOs(MOEnergy, CMatrixTTemp, CMatrixT, NumberOfOrbitals);
    CMatrix = Sminus * CMatrixT;
    //Normalization of C matrix
    NormalizationC(CMatrix, OverlapIntegralMatrix);
    //Calculate the density matrix
    Matrix<double, Dynamic, Dynamic> P;
    P.resize(NumberOfOrbitals, NumberOfOrbitals);
    for (int i = 0; i < NumberOfOrbitals; i++)
    {
        for (int j = 0; j < NumberOfOrbitals; j++)
        {
            double Ptu = 0;
            for (int a = 0; a < Electrons / 2; a++)
            {
                Ptu += CMatrix(i, a) * CMatrix(j, a);
            }
            Ptu *= 2.0;
            P(i, j) = Ptu;
        }
    }
    MOEnergy.clear();
    HFout << "******** Overlap Matrix ********" << endl;
    HFout << OverlapIntegralMatrix << endl;
    HFout << endl;
    HFout << "******** Kinetic Matrix ********" << endl;
    HFout << KineticEnergyMatrix << endl;
    HFout << endl;
    HFout << "******** Potential Energy Matrix ********" << endl;
    HFout << NuclearAttractionMatrix << endl;
    HFout << endl;
    HFout << "******** HCore Matrix ********" << endl;
    HFout << HCore << endl;
    HFout << endl;
    HFout << "Nuclear repulsive energy = " << RepulsiveEnergy << endl;
    //SCF starts
    int CountSCFCycle = 1;
    vector<double> EnergyHistory;
    int IfConverge = 0;
    double EnergyDifference = 100;
    for (int i = 0; i < 3; i++)
    {
        HFout << endl;
    }
    while (IfConverge == 0)
    {
        HFout << "SCF cycle = " << CountSCFCycle << endl;
        HFout << endl;
        HFout << "******** Density Matrix ********" << endl;
        HFout << P << endl;
        HFout << endl;
        //Calculate the Fock matrix
        for (int i = 0; i < NumberOfOrbitals; i++)
        {
            for (int j = 0; j < NumberOfOrbitals; j++)
            {
                double KandJIntegral = 0;
                for (int tt = 0; tt < NumberOfOrbitals; tt++)
                {
                    for (int uu = 0; uu < NumberOfOrbitals; uu++)
                    {
                        KandJIntegral += P(tt, uu) * (gabcd[i][j][tt][uu] - 0.5 * gabcd[i][uu][tt][j]);
                    }
                }
                Fock(i, j) = HCore(i, j) + KandJIntegral;
                if (fabs(Fock(i, j)) < pow(10, -9))
                {
                    Fock(i, j) = 0;
                }
            }
        }
        HFout << "******** Fock Matrix ********" << endl;
        HFout << Fock << endl;
        HFout << endl;
        double ElectronEnergy = 0;
        for (int i = 0; i < NumberOfOrbitals; i++)
        {
            for (int j = 0; j < NumberOfOrbitals; j++)
            {
                ElectronEnergy += 0.5 * P(i, j) * (Fock(i, j) + HCore(i, j));
            }
        }
        Energy = ElectronEnergy + RepulsiveEnergy;
        HFout << "The Energy of molecule is: " << Energy << endl;
        HFout << endl;
        if (CountSCFCycle > 1)
        {
            EnergyDifference = Energy - EnergyHistory[EnergyHistory.size() - 1];
        }
        EnergyHistory.push_back(Energy);
        //Calculate the FockT matrix (F')
        FockT = Sminus * Fock * Sminus;
        //Diagonalize the FockT matrix
        EigenSolver<MatrixXd> EnergySolution(FockT);
        CMatrixTTemp = EnergySolution.eigenvectors().real();
        for (int i = 0; i < EnergySolution.eigenvalues().size(); i++)
        {
            pair<int, double> tp(i, EnergySolution.eigenvalues()[i].real());
            MOEnergy.push_back(tp);
        }
        stable_sort(MOEnergy.begin(), MOEnergy.end(), CompareEnergy);
        ArrangeMOs(MOEnergy, CMatrixTTemp, CMatrixT, NumberOfOrbitals);
        //Update the C Matrix and density matrix
        CMatrix = Sminus * CMatrixT;
        NormalizationC(CMatrix, OverlapIntegralMatrix);
        HFout << "******** Molecular Orbitals ********" << endl;
        HFout << endl;
        HFout << "The Orbitals are ranged in the sequence below:" << endl;
        for (int k = 0; k < NumberOfOrbitals; k++)
        {
            string AOName;
            PrintOrbital(OrbitalSequence[k].SerialNumber, OrbitalSequence[k].n, OrbitalSequence[k].l, OrbitalSequence[k].m, AOName, inp);
            if (k < NumberOfOrbitals - 1)
            {
                HFout << AOName << " ";
            }
            else
            {
                HFout << AOName << endl;
            }
        }
        HFout << endl;
        HFout << "The energy of molecular orbitals are:" << endl;
        for (int i = 0; i < NumberOfOrbitals; i++)
        {
            HFout << "MO " << i + 1 << ": " << MOEnergy[i].second << " ";
            HFout << "Number of electrons: ";
            if (i < Electrons / 2)
            {
                HFout << 2 << endl;
            }
            else
            {
                HFout << 0 << endl;
            }
            for (int j = 0; j < NumberOfOrbitals; j++)
            {
                double AOCo = CMatrix(j, i);
                if (fabs(AOCo) < pow(10, -5))
                {
                    AOCo = 0;
                }
                if (j == 0)
                {
                    HFout << AOCo;
                }
                else
                {
                    HFout << " " << AOCo;
                }
                if (j == NumberOfOrbitals - 1)
                {
                    HFout << endl;
                }
            }
        }
        for (int i = 0; i < 3; i++)
        {
            HFout << endl;
        }
        for (int i = 0; i < NumberOfOrbitals; i++)
        {
            for (int j = 0; j < NumberOfOrbitals; j++)
            {
                double Ptu = 0;
                for (int a = 0; a < Electrons / 2; a++)
                {
                    Ptu += CMatrix(i, a) * CMatrix(j, a);
                }
                Ptu *= 2;
                P(i, j) = Ptu;
            }
        }
        if (fabs(EnergyDifference) < pow(10, -9) || CountSCFCycle > 30)
        {
            IfConverge = 1;
        }
        MOEnergy.clear();
        CountSCFCycle++;
    }
}