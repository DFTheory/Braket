#include "Controller.hpp"
using namespace std;

double HermiteCo(int i, int j, int t, double XPA, double XPB, double p, double KAB)
{
    double result = 0;
    if (t < 0 || t > i + j)
    {
        result = 0;
    }
    else
    {
        if (i == 0 && j == 0)
        {
            result = KAB;
        }
        else
        {
            if (i > 0)
            {
                if (j == 0)
                {
                    result = HermiteCo(i - 1, j, t - 1, XPA, XPB, p, KAB) / (2 * p) + XPA * HermiteCo(i - 1, j, t, XPA, XPB, p, KAB) + (t + 1) * HermiteCo(i - 1, j, t + 1, XPA, XPB, p, KAB);
                }
                else
                {
                    result = HermiteCo(i, j - 1, t - 1, XPA, XPB, p, KAB) / (2 * p) + XPB * HermiteCo(i, j - 1, t, XPA, XPB, p, KAB) + (t + 1) * HermiteCo(i, j - 1, t + 1, XPA, XPB, p, KAB);
                }
            }
            else if (i == 0)
            {
                result = HermiteCo(i, j - 1, t - 1, XPA, XPB, p, KAB) / (2 * p) + XPB * HermiteCo(i, j - 1, t, XPA, XPB, p, KAB) + (t + 1) * HermiteCo(i, j - 1, t + 1, XPA, XPB, p, KAB);
            }
        }
    }
    return result;
}

double SIntegral(double rA, double rB, double aA, double aB, int i, int j) //This function is only in one dimension and should be expanded to 3 dimension.
{
    double p, KAB, mu, AB, PA, PB;
    p = aA + aB;
    AB = rA - rB;
    PA = (aA * rA + aB * rB) / p - rA;
    PB = (aA * rA + aB * rB) / p - rB;
    mu = aA * aB / (aA + aB);
    KAB = exp(-mu * AB * AB);
    return HermiteCo(i, j, 0, PA, PB, p, KAB) * sqrt(M_PI / p);
}

double OrbitalSIntegral(Atom p, Atom q, Orbital A, Orbital B)
{
    double result = 0;
    for (int i = 0; i < A.Expression.size(); i++)
    {
        for (int j = 0; j < B.Expression.size(); j++)
        {
            double xIntegral, yIntegral, zIntegral;
            xIntegral = SIntegral(p.Cartesian[0], q.Cartesian[0], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].xi, B.Expression[j].xi);
            yIntegral = SIntegral(p.Cartesian[1], q.Cartesian[1], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].yj, B.Expression[j].yj);
            zIntegral = SIntegral(p.Cartesian[2], q.Cartesian[2], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].zk, B.Expression[j].zk);
            result += xIntegral * yIntegral * zIntegral * A.Expression[i].ecofficient * B.Expression[j].ecofficient;
        }
    }
    return result;
}

double OrbitalKineticEnergyIntegral(Atom p, Atom q, Orbital A, Orbital B)
{
    double result = 0;
    for (int i = 0; i < A.Expression.size(); i++)
    {
        for (int j = 0; j < B.Expression.size(); j++)
        {
            double xIntegral, yIntegral, zIntegral;
            double xSIntegral, ySIntegral, zSIntegral;
            xSIntegral = SIntegral(p.Cartesian[0], q.Cartesian[0], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].xi, B.Expression[j].xi);
            ySIntegral = SIntegral(p.Cartesian[1], q.Cartesian[1], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].yj, B.Expression[j].yj);
            zSIntegral = SIntegral(p.Cartesian[2], q.Cartesian[2], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].zk, B.Expression[j].zk);
            xIntegral = 4 * pow(B.Expression[j].ealpha, 2) * SIntegral(p.Cartesian[0], q.Cartesian[0], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].xi, B.Expression[j].xi + 2) - 2 * B.Expression[j].ealpha * (2 * B.Expression[j].xi + 1) * xSIntegral + B.Expression[j].xi * (B.Expression[j].xi - 1) * SIntegral(p.Cartesian[0], q.Cartesian[0], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].xi, B.Expression[j].xi - 2);
            yIntegral = 4 * pow(B.Expression[j].ealpha, 2) * SIntegral(p.Cartesian[1], q.Cartesian[1], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].yj, B.Expression[j].yj + 2) - 2 * B.Expression[j].ealpha * (2 * B.Expression[j].yj + 1) * ySIntegral + B.Expression[j].yj * (B.Expression[j].yj - 1) * SIntegral(p.Cartesian[1], q.Cartesian[1], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].yj, B.Expression[j].yj - 2);
            zIntegral = 4 * pow(B.Expression[j].ealpha, 2) * SIntegral(p.Cartesian[2], q.Cartesian[2], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].zk, B.Expression[j].zk + 2) - 2 * B.Expression[j].ealpha * (2 * B.Expression[j].zk + 1) * zSIntegral + B.Expression[j].zk * (B.Expression[j].zk - 1) * SIntegral(p.Cartesian[2], q.Cartesian[2], A.Expression[i].ealpha, B.Expression[j].ealpha, A.Expression[i].zk, B.Expression[j].zk - 2);
            result += (A.Expression[i].ecofficient * B.Expression[j].ecofficient) * (xIntegral * ySIntegral * zSIntegral + xSIntegral * yIntegral * zSIntegral + xSIntegral * ySIntegral * zIntegral);
        }
    }
    return result * (-0.5);
}

int DoubleFactorial(int x)
{
    if (x == -1 || x == 0)
    {
        return 1;
    }
    else
    {
        int result = 1;
        if (x % 2 == 0)
        {
            for (int i = 2; i <= x; i += 2)
            {
                result *= i;
            }
        }
        else
        {
            for (int i = 1; i <= x; i += 2)
            {
                result *= i;
            }
        }
        return result;
    }
}

double BoysFunction(int n, double x)
{
    double result = 0;
    if (x > 35)
    {
        result = (DoubleFactorial(2 * n - 1) / (pow(2, n + 1))) * sqrt(M_PI / (pow(x, 2 * n + 1)));
    }
    else
    {
        double partA = exp(-x);
        double partD = 0;
        for (int i = 0; i < 100; i++)
        {
            double partB = pow(2 * x, i);
            double partC = 1;
            for (int j = 0; j <= i; j++)
            {
                partC *= (2 * n + 2 * j + 1);
            }
            partD += partB / partC;
        }
        result = partA * partD;
    }
    return result;
}

double BoysFunctionRecurrence(int n, double F, double x)
{
    double result = (2 * x * F + exp(-x)) / (2 * n - 1);
    return result;
}

double Rtuv(int t, int u, int v, int n, double p, double Ax, double Ay, double Az, double FVal[])
{
    double result;
    if (t < 0 || u < 0 || v < 0 || n < 0)
    {
        result = 0;
    }
    else
    {
        if (t == 0 && u == 0 && v == 0)
        {
            result = pow(-2 * p, n) * FVal[n];
        }
        else
        {
            if (v > 0)
            {
                result = (v - 1) * Rtuv(t, u, v - 2, n + 1, p, Ax, Ay, Az, FVal) + Az * Rtuv(t, u, v - 1, n + 1, p, Ax, Ay, Az, FVal);
            }
            else if (v == 0)
            {
                if (u > 0)
                {
                    result = (u - 1) * Rtuv(t, u - 2, v, n + 1, p, Ax, Ay, Az, FVal) + Ay * Rtuv(t, u - 1, v, n + 1, p, Ax, Ay, Az, FVal);
                }
                else if (u == 0)
                {
                    if (t > 0)
                    {
                        result = (t - 1) * Rtuv(t - 2, u, v, n + 1, p, Ax, Ay, Az, FVal) + Ax * Rtuv(t - 1, u, v, n + 1, p, Ax, Ay, Az, FVal);
                    }
                }
            }
        }
    }
    return result;
}

double NuclearAttractionIntegral(Atom a, Atom b, GaussTypeFunction x, GaussTypeFunction y, vector<Atom> inp)
{
    double p = x.ealpha + y.ealpha;
    double result = 0;
    double px, py, pz;
    px = (x.ealpha * a.Cartesian[0] + y.ealpha * b.Cartesian[0]) / p;
    py = (x.ealpha * a.Cartesian[1] + y.ealpha * b.Cartesian[1]) / p;
    pz = (x.ealpha * a.Cartesian[2] + y.ealpha * b.Cartesian[2]) / p;
    double mu = x.ealpha * y.ealpha / p;
    double XAB = a.Cartesian[0] - b.Cartesian[0];
    double YAB = a.Cartesian[1] - b.Cartesian[1];
    double ZAB = a.Cartesian[2] - b.Cartesian[2];
    double KABx = exp(-mu * XAB * XAB);
    double KABy = exp(-mu * YAB * YAB);
    double KABz = exp(-mu * ZAB * ZAB);
    int total_t = x.xi + y.xi;
    int total_u = x.yj + y.yj;
    int total_v = x.zk + y.zk;
    int total = total_t + total_u + total_v;
    for (int i = 0; i < inp.size(); i++)
    {
        double temp = 0;
        double RCP;
        double Cx, Cy, Cz;
        Cx = inp[i].Cartesian[0], Cy = inp[i].Cartesian[1], Cz = inp[i].Cartesian[2];
        RCP = sqrt(pow(Cx - px, 2) + pow(Cy - py, 2) + pow(Cz - pz, 2));
        double XPC = px - Cx, YPC = py - Cy, ZPC = pz - Cz;
        double F = BoysFunction(total, p * RCP * RCP);
        double FValues[1000];
        FValues[total] = F;
        for (int j = total - 1; j >= 0; j--)
        {
            FValues[j] = BoysFunctionRecurrence(j + 1, FValues[j + 1], p * RCP * RCP);
        }
        for (int t = 0; t <= total_t; t++)
        {
            for (int u = 0; u <= total_u; u++)
            {
                for (int v = 0; v <= total_v; v++)
                {
                    double Eijt = HermiteCo(x.xi, y.xi, t, px - a.Cartesian[0], px - b.Cartesian[0], p, KABx);
                    double Eklu = HermiteCo(x.yj, y.yj, u, py - a.Cartesian[1], py - b.Cartesian[1], p, KABy);
                    double Emnv = HermiteCo(x.zk, y.zk, v, pz - a.Cartesian[2], pz - b.Cartesian[2], p, KABz);
                    temp += Eijt * Eklu * Emnv * Rtuv(t, u, v, 0, p, XPC, YPC, ZPC, FValues);
                }
            }
        }
        temp *= -inp[i].AtomSerialNumber;
        result += temp;
    }
    result *= (2 * M_PI / p);
    return result;
}

double OrbitalNuclearAttractionIntegral(Atom a, Atom b, Orbital A, Orbital B, vector<Atom> inpu)
{
    double result = 0;
    for (int i = 0; i < A.Expression.size(); i++)
    {
        for (int j = 0; j < B.Expression.size(); j++)
        {
            result += NuclearAttractionIntegral(a, b, A.Expression[i], B.Expression[j], inpu) * A.Expression[i].ecofficient * B.Expression[j].ecofficient;
        }
    }
    return result;
}

double FourCenterTwoElectronIntegral(Atom A, Atom B, Atom C, Atom D, GaussTypeFunction a, GaussTypeFunction b, GaussTypeFunction c, GaussTypeFunction d)
{
    double result = 0;
    double p = a.ealpha + b.ealpha;
    double q = c.ealpha + d.ealpha;
    double Alpha = p * q / (p + q);
    double TheCo = (2 * pow(M_PI, 2.5)) / (p * q * sqrt(p + q));
    double FValues[1000];
    double RPQ, Px, Py, Pz, Qx, Qy, Qz, PQx, PQy, PQz;
    Px = (a.ealpha * A.Cartesian[0] + b.ealpha * B.Cartesian[0]) / p;
    Py = (a.ealpha * A.Cartesian[1] + b.ealpha * B.Cartesian[1]) / p;
    Pz = (a.ealpha * A.Cartesian[2] + b.ealpha * B.Cartesian[2]) / p;
    Qx = (c.ealpha * C.Cartesian[0] + d.ealpha * D.Cartesian[0]) / q;
    Qy = (c.ealpha * C.Cartesian[1] + d.ealpha * D.Cartesian[1]) / q;
    Qz = (c.ealpha * C.Cartesian[2] + d.ealpha * D.Cartesian[2]) / q;
    PQx = Px - Qx;
    PQy = Py - Qy;
    PQz = Pz - Qz;
    RPQ = sqrt(PQx * PQx + PQy * PQy + PQz * PQz);
    double mup, muq, KABx, KABy, KABz, KCDx, KCDy, KCDz;
    mup = a.ealpha * b.ealpha / p;
    muq = c.ealpha * d.ealpha / q;
    KABx = exp(-mup * pow(A.Cartesian[0] - B.Cartesian[0], 2));
    KABy = exp(-mup * pow(A.Cartesian[1] - B.Cartesian[1], 2));
    KABz = exp(-mup * pow(A.Cartesian[2] - B.Cartesian[2], 2));
    KCDx = exp(-muq * pow(C.Cartesian[0] - D.Cartesian[0], 2));
    KCDy = exp(-muq * pow(C.Cartesian[1] - D.Cartesian[1], 2));
    KCDz = exp(-muq * pow(C.Cartesian[2] - D.Cartesian[2], 2));
    int total_tI, total_uI, total_vI, total_tII, total_uII, total_vII, total;
    total_tI = a.xi + b.xi;
    total_uI = a.yj + b.yj;
    total_vI = a.zk + b.zk;
    total_tII = c.xi + d.xi;
    total_uII = c.yj + d.yj;
    total_vII = c.zk + d.zk;
    total = total_tI + total_tII + total_uI + total_uII + total_vI + total_vII;
    FValues[total] = BoysFunction(total, Alpha * RPQ * RPQ);
    for (int i = total - 1; i >= 0; i--)
    {
        FValues[i] = BoysFunctionRecurrence(i + 1, FValues[i + 1], Alpha * RPQ * RPQ);
    }
    for (int tI = 0; tI <= total_tI; tI++)
    {
        for (int uI = 0; uI <= total_uI; uI++)
        {
            for (int vI = 0; vI <= total_vI; vI++)
            {
                double EtI, EuI, EvI;
                EtI = HermiteCo(a.xi, b.xi, tI, Px - A.Cartesian[0], Px - B.Cartesian[0], p, KABx);
                EuI = HermiteCo(a.yj, b.yj, uI, Py - A.Cartesian[1], Py - B.Cartesian[1], p, KABy);
                EvI = HermiteCo(a.zk, b.zk, vI, Pz - A.Cartesian[2], Pz - B.Cartesian[2], p, KABz);
                double EI = EtI * EuI * EvI;
                for (int tII = 0; tII <= total_tII; tII++)
                {
                    for (int uII = 0; uII <= total_uII; uII++)
                    {
                        for (int vII = 0; vII <= total_vII; vII++)
                        {
                            double EtII, EuII, EvII;
                            EtII = HermiteCo(c.xi, d.xi, tII, Qx - C.Cartesian[0], Qx - D.Cartesian[0], q, KCDx);
                            EuII = HermiteCo(c.yj, d.yj, uII, Qy - C.Cartesian[1], Qy - D.Cartesian[1], q, KCDy);
                            EvII = HermiteCo(c.zk, d.zk, vII, Qz - C.Cartesian[2], Qz - D.Cartesian[2], q, KCDz);
                            double EII = EtII * EuII * EvII;
                            result += EI * EII * Rtuv(tI + tII, uI + uII, vI + vII, 0, Alpha, PQx, PQy, PQz, FValues) * pow(-1, tII + uII + vII);
                        }
                    }
                }
            }
        }
    }
    result *= TheCo;
    return result;
}

double OrbitalFCTEIntegral(Atom A, Atom B, Atom C, Atom D, Orbital a, Orbital b, Orbital c, Orbital d) //FCTE means four center two electrons.
{
    double result = 0;
    for (int i = 0; i < a.Expression.size(); i++)
    {
        for (int j = 0; j < b.Expression.size(); j++)
        {
            for (int p = 0; p < c.Expression.size(); p++)
            {
                for (int q = 0; q < d.Expression.size(); q++)
                {
                    double coef = a.Expression[i].ecofficient * b.Expression[j].ecofficient * c.Expression[p].ecofficient * d.Expression[q].ecofficient;
                    result += coef * FourCenterTwoElectronIntegral(A, B, C, D, a.Expression[i], b.Expression[j], c.Expression[p], d.Expression[q]);
                }
            }
        }
    }
    return result;
}