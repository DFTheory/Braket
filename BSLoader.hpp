#include "OrbitalBank.hpp"
using namespace std;

int AtomSelector(string atom) //return the serial number of an element
{
    int findit = Elements.count(atom);
    if (findit == 0)
    {
        cout << "Atom not found or defined" << endl;
        abort();
    }
    else
    {
        return Elements[atom];
    }
}

double GsStyleToDouble(string a)
{
    string part1, part2;
    double p1, p2, re;
    part1 = a.substr(0, a.length() - 4);
    part2 = a.substr(a.length() - 3);
    p1 = atof(part1.c_str());
    p2 = atof(part2.c_str());
    return p1 * pow(10, p2);
}

int FunctionNumber(int Serial, int bs)
{
    int FunctionNumbers;
    if (bs == 1)
    {
        if (Serial < 2)
        {
            FunctionNumbers = 1;
        }
        else if (Serial >= 2 && Serial < 10)
        {
            FunctionNumbers = 2;
        }
        else if (Serial >= 10 && Serial < 18)
        {
            FunctionNumbers = 3;
        }
        else if (Serial >= 18 && Serial < 20)
        {
            FunctionNumbers = 4;
        }
        else if (Serial >= 20 && Serial < 36)
        {
            FunctionNumbers = 5;
        }
        else if (Serial >= 36 && Serial < 38)
        {
            FunctionNumbers = 6;
        }
        else if (Serial >= 38 && Serial < 54)
        {
            FunctionNumbers = 7;
        }
    }
    return FunctionNumbers;
}

void BuildBSBank(vector<BSUnit> bsb, int bs)
{
    if (bs == 1)
    {
        fstream iptsto3g("sto-3g.1.gbs", ios::in);
        int TotalAtoms;
        iptsto3g >> TotalAtoms;
        for (int i = 0; i < TotalAtoms; i++)
        {
            BSUnit TEMP;
            string atom, endline;
            int a;
            iptsto3g >> atom >> a;
            int BasisFunctionNumbers = FunctionNumber(i, 1); //Define the number of orbital sets
            for (int j = 0; j < BasisFunctionNumbers; j++)
            {
                string TypeOfOrbital;
                int Amount;
                double ScaleFactor;
                iptsto3g >> TypeOfOrbital >> Amount >> ScaleFactor;
                vector<vector<double>> tempS, tempP, tempD, tempF;
                for (int k = 0; k < Amount; k++)
                {
                    vector<double> tp;
                    if (TypeOfOrbital.compare("SP") == 0)
                    {
                        string a1, c1, c2;
                        iptsto3g >> a1 >> c1 >> c2;
                        double a1d = GsStyleToDouble(a1), csd = GsStyleToDouble(c1), cpd = GsStyleToDouble(c2);
                        tp.push_back(a1d), tp.push_back(csd);
                        tempS.push_back(tp);
                        tp.pop_back(), tp.push_back(cpd);
                        tempP.push_back(tp);
                    }
                    else
                    {
                        string alpha, co;
                        iptsto3g >> alpha >> co;
                        double alphad = GsStyleToDouble(alpha), cod = GsStyleToDouble(co);
                        tp.push_back(alphad), tp.push_back(cod);
                        if (TypeOfOrbital.compare("S") == 0)
                        {
                            tempS.push_back(tp);
                        }
                        else if (TypeOfOrbital.compare("P") == 0)
                        {
                            tempP.push_back(tp);
                        }
                        else if (TypeOfOrbital.compare("D") == 0)
                        {
                            tempD.push_back(tp);
                        }
                        else if (TypeOfOrbital.compare("F") == 0)
                        {
                            tempF.push_back(tp);
                        }
                    }
                }
                if (tempS.size() > 0)
                {
                    TEMP.S.push_back(tempS);
                }
                if (tempP.size() > 0)
                {
                    TEMP.P.push_back(tempP);
                }
                if (tempD.size() > 0)
                {
                    TEMP.D.push_back(tempD);
                }
                if (tempF.size() > 0)
                {
                    TEMP.F.push_back(tempF);
                }
            }
            iptsto3g >> endline;
            BSBank.push_back(TEMP);
        }
    }
}

void OrbitalMaker(int bs, int ProtonNumber, Atom &example) //This function helps connect each orbitals with basis functions.
{
    if (bs == 1)
    {
        ProtonNumber -= 1;
        vector<vector<Orbital>> shell;
        vector<Orbital> tp;
        //n = 1 aka as 1s orbital.
        Orbital OneSobt;
        for (int i = 0; i < BSBank[ProtonNumber].S[0].size(); i++)
        {
            GaussTypeFunction gtf;
            gtf.xi = 0, gtf.yj = 0, gtf.zk = 0;
            gtf.ealpha = BSBank[ProtonNumber].S[0][i][0];
            gtf.ecofficient = BSBank[ProtonNumber].S[0][i][1] * NmlCo(0, 0, 0, gtf.ealpha);
            OneSobt.Expression.push_back(gtf);
        }
        tp.push_back(OneSobt);
        shell.push_back(tp);
        example.AtomOrbitals.push_back(shell);
        tp.clear();
        shell.clear();
        //n = 2, this shell includes one 2s orbital and three 2p orbitals.
        //2s first
        if (BSBank[ProtonNumber].S.size() > 1)
        {
            Orbital TwoSobt;
            for (int i = 0; i < BSBank[ProtonNumber].S[1].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.xi = 0, gtf.yj = 0, gtf.zk = 0;
                gtf.ealpha = BSBank[ProtonNumber].S[1][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].S[1][i][1] * NmlCo(0, 0, 0, gtf.ealpha);
                TwoSobt.Expression.push_back(gtf);
            }
            tp.push_back(TwoSobt);
            shell.push_back(tp);
            tp.clear();
        }
        //2p
        if (BSBank[ProtonNumber].P.size() > 0)
        {
            Orbital TwoPx, TwoPy, TwoPz;
            for (int i = 0; i < BSBank[ProtonNumber].P[0].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.ealpha = BSBank[ProtonNumber].P[0][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].P[0][i][1] * NmlCo(1, 0, 0, gtf.ealpha);
                gtf.xi = 1, gtf.yj = 0, gtf.zk = 0;
                TwoPx.Expression.push_back(gtf);
                gtf.xi = 0, gtf.yj = 1;
                TwoPy.Expression.push_back(gtf);
                gtf.yj = 0, gtf.zk = 1;
                TwoPz.Expression.push_back(gtf);
            }
            tp.push_back(TwoPx), tp.push_back(TwoPy), tp.push_back(TwoPz);
            shell.push_back(tp);
            tp.clear();
        }
        if (shell.size() > 0)
        {
            example.AtomOrbitals.push_back(shell);
            shell.clear();
        }
        //n = 3
        //3s
        if (BSBank[ProtonNumber].S.size() > 2)
        {
            Orbital ThreeSobt;
            for (int i = 0; i < BSBank[ProtonNumber].S[2].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.xi = 0, gtf.yj = 0, gtf.zk = 0;
                gtf.ealpha = BSBank[ProtonNumber].S[2][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].S[2][i][1] * NmlCo(0, 0, 0, gtf.ealpha);
                ThreeSobt.Expression.push_back(gtf);
            }
            tp.push_back(ThreeSobt);
            shell.push_back(tp);
            tp.clear();
        }
        //3p
        if (BSBank[ProtonNumber].P.size() > 1)
        {
            Orbital ThreePx, ThreePy, ThreePz;
            for (int i = 0; i < BSBank[ProtonNumber].P[1].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.ealpha = BSBank[ProtonNumber].P[1][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].P[1][i][1] * NmlCo(1, 0, 0, gtf.ealpha);
                gtf.xi = 1, gtf.yj = 0, gtf.zk = 0;
                ThreePx.Expression.push_back(gtf);
                gtf.xi = 0, gtf.yj = 1;
                ThreePy.Expression.push_back(gtf);
                gtf.yj = 0, gtf.zk = 1;
                ThreePz.Expression.push_back(gtf);
            }
            tp.push_back(ThreePx), tp.push_back(ThreePy), tp.push_back(ThreePz);
            shell.push_back(tp);
            tp.clear();
        }
        //3d
        if (BSBank[ProtonNumber].D.size() > 0)
        {
            Orbital ThreeDzero, ThreeDpOne, ThreeDmOne, ThreeDmTwo, ThreeDpTwo;
            for (int i = 0; i < BSBank[ProtonNumber].D[0].size(); i++)
            {
                GaussTypeFunction gtfa, gtfb, gtfc;
                gtfa.ealpha = BSBank[ProtonNumber].D[0][i][0];
                gtfb.ealpha = gtfa.ealpha, gtfc.ealpha = gtfa.ealpha;
                double c = BSBank[ProtonNumber].D[0][i][1];
                gtfa.ecofficient = c * (-0.5) * NmlCo(2, 0, 0, gtfa.ealpha);
                gtfa.xi = 2, gtfa.yj = 0, gtfa.zk = 0;
                gtfb.ecofficient = c * (-0.5) * NmlCo(0, 2, 0, gtfb.ealpha);
                gtfb.xi = 0, gtfb.yj = 2, gtfb.zk = 0;
                gtfc.ecofficient = c * NmlCo(0, 0, 2, gtfc.ealpha);
                gtfc.xi = 0, gtfc.yj = 0, gtfc.zk = 2;
                ThreeDzero.Expression.push_back(gtfa), ThreeDzero.Expression.push_back(gtfb), ThreeDzero.Expression.push_back(gtfc);
                gtfa.ecofficient = c * NmlCo(1, 0, 1, gtfa.ealpha);
                gtfa.xi = 1, gtfa.yj = 0, gtfa.zk = 1;
                ThreeDpOne.Expression.push_back(gtfa);
                gtfa.xi = 0, gtfa.yj = 1;
                ThreeDmOne.Expression.push_back(gtfa);
                gtfa.xi = 1, gtfa.zk = 0;
                ThreeDmTwo.Expression.push_back(gtfa);
                gtfa.ecofficient = c * (sqrt(3) / 2) * NmlCo(2, 0, 0, gtfa.ealpha);
                gtfb.ecofficient = -gtfa.ecofficient;
                gtfa.xi = 2, gtfa.yj = 0, gtfa.zk = 0;
                gtfb.xi = 0, gtfb.yj = 2, gtfb.zk = 0;
                ThreeDpTwo.Expression.push_back(gtfa), ThreeDpTwo.Expression.push_back(gtfb);
            }
            tp.push_back(ThreeDzero), tp.push_back(ThreeDpOne), tp.push_back(ThreeDmOne), tp.push_back(ThreeDpTwo), tp.push_back(ThreeDmTwo);
            shell.push_back(tp);
            tp.clear();
        }
        if (shell.size() > 0)
        {
            example.AtomOrbitals.push_back(shell);
            shell.clear();
        }
        //n = 4
        //4s
        if (BSBank[ProtonNumber].S.size() > 3)
        {
            Orbital FourS;
            for (int i = 0; i < BSBank[ProtonNumber].S[3].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.ealpha = BSBank[ProtonNumber].S[3][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].S[3][i][1] * NmlCo(0, 0, 0, gtf.ealpha);
                gtf.xi = 0, gtf.yj = 0, gtf.zk = 0;
                FourS.Expression.push_back(gtf);
            }
            tp.push_back(FourS);
            shell.push_back(tp);
            tp.clear();
        }
        //4p
        if (BSBank[ProtonNumber].P.size() > 2)
        {
            Orbital FourPx, FourPy, FourPz;
            for (int i = 0; i < BSBank[ProtonNumber].P[2].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.ealpha = BSBank[ProtonNumber].P[2][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].P[2][i][1] * NmlCo(1, 0, 0, gtf.ealpha);
                gtf.xi = 1, gtf.yj = 0, gtf.zk = 0;
                FourPx.Expression.push_back(gtf);
                gtf.xi = 0, gtf.yj = 1;
                FourPy.Expression.push_back(gtf);
                gtf.yj = 0, gtf.zk = 1;
                FourPz.Expression.push_back(gtf);
            }
            tp.push_back(FourPx), tp.push_back(FourPy), tp.push_back(FourPz);
            shell.push_back(tp);
            tp.clear();
        }
        //4d
        if (BSBank[ProtonNumber].D.size() > 1)
        {
            Orbital FourDzero, FourDpOne, FourDmOne, FourDmTwo, FourDpTwo;
            for (int i = 0; i < BSBank[ProtonNumber].D[1].size(); i++)
            {
                GaussTypeFunction gtfa, gtfb, gtfc;
                gtfa.ealpha = BSBank[ProtonNumber].D[1][i][0];
                gtfb.ealpha = gtfa.ealpha, gtfc.ealpha = gtfa.ealpha;
                double c = BSBank[ProtonNumber].D[1][i][1];
                gtfa.ecofficient = c * (-0.5) * NmlCo(2, 0, 0, gtfa.ealpha);
                gtfa.xi = 2, gtfa.yj = 0, gtfa.zk = 0;
                gtfb.ecofficient = c * (-0.5) * NmlCo(0, 2, 0, gtfb.ealpha);
                gtfb.xi = 0, gtfb.yj = 2, gtfb.zk = 0;
                gtfc.ecofficient = c * NmlCo(0, 0, 2, gtfc.ealpha);
                gtfc.xi = 0, gtfc.yj = 0, gtfc.zk = 2;
                FourDzero.Expression.push_back(gtfa), FourDzero.Expression.push_back(gtfb), FourDzero.Expression.push_back(gtfc);
                gtfa.ecofficient = c * NmlCo(1, 0, 1, gtfa.ealpha);
                gtfa.xi = 1, gtfa.yj = 0, gtfa.zk = 1;
                FourDpOne.Expression.push_back(gtfa);
                gtfa.xi = 0, gtfa.yj = 1;
                FourDmOne.Expression.push_back(gtfa);
                gtfa.xi = 1, gtfa.zk = 0;
                FourDmTwo.Expression.push_back(gtfa);
                gtfa.ecofficient = c * (sqrt(3) / 2) * NmlCo(2, 0, 0, gtfa.ealpha);
                gtfb.ecofficient = -gtfa.ecofficient;
                gtfa.xi = 2, gtfa.yj = 0, gtfa.zk = 0;
                gtfb.xi = 0, gtfb.yj = 2, gtfb.zk = 0;
                FourDpTwo.Expression.push_back(gtfa), FourDpTwo.Expression.push_back(gtfb);
            }
            tp.push_back(FourDzero), tp.push_back(FourDpOne), tp.push_back(FourDmOne), tp.push_back(FourDpTwo), tp.push_back(FourDmTwo);
            shell.push_back(tp);
            tp.clear();
        }
        if (shell.size() > 0)
        {
            example.AtomOrbitals.push_back(shell);
            shell.clear();
        }
        //4f orbitals are not considered in STO-3g basis sets
        //n = 5
        //5s
        if (BSBank[ProtonNumber].S.size() > 4)
        {
            Orbital FiveS;
            for (int i = 0; i < BSBank[ProtonNumber].S[4].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.xi = 0, gtf.yj = 0, gtf.zk = 0;
                gtf.ealpha = BSBank[ProtonNumber].S[4][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].S[4][i][1] * NmlCo(0, 0, 0, gtf.ealpha);
                FiveS.Expression.push_back(gtf);
            }
            tp.push_back(FiveS);
            shell.push_back(tp);
            tp.clear();
        }
        //5p
        if (BSBank[ProtonNumber].P.size() > 3)
        {
            Orbital FivePx, FivePy, FivePz;
            for (int i = 0; i < BSBank[ProtonNumber].P[3].size(); i++)
            {
                GaussTypeFunction gtf;
                gtf.ealpha = BSBank[ProtonNumber].P[3][i][0];
                gtf.ecofficient = BSBank[ProtonNumber].P[3][i][1] * NmlCo(1, 0, 0, gtf.ealpha);
                gtf.xi = 1, gtf.yj = 0, gtf.zk = 0;
                FivePx.Expression.push_back(gtf);
                gtf.xi = 0, gtf.yj = 1;
                FivePy.Expression.push_back(gtf);
                gtf.yj = 0, gtf.zk = 1;
                FivePz.Expression.push_back(gtf);
            }
            tp.push_back(FivePx), tp.push_back(FivePy), tp.push_back(FivePz);
            shell.push_back(tp);
            tp.clear();
        }
        if (shell.size() > 0)
        {
            example.AtomOrbitals.push_back(shell);
            shell.clear();
        }
        //To be continued...
    }
}