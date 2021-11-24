#include "BSLoader.hpp"
using namespace std;

vector<vector<vector<double>>> BasisSetsBank;
int JobControl(string job)
{
    string SPE("SPE");
    string opt("OPT");
    transform(job.begin(), job.end(), job.begin(), ::toupper);
    if (job.compare(SPE) == 0)
    {
        return 1;
    }
    else if (job.compare(opt) == 0)
    {
        return 2;
    }
    else
    {
        cout << "Job type not developed yet or wrong" << endl;
        abort();
    }
}

int BasisSetsControl(string bs)
{
    string STO3G("STO-3G");
    transform(bs.begin(), bs.end(), bs.begin(), ::toupper);
    if (bs.compare(STO3G) == 0)
    {
        return 1;
    }
    else
    {
        cout << "Basis Sets not included or defined" << endl;
        abort();
    }
}

int MethodControl(string mtd)
{
    string HF("HF");
    string B3LYP("B3LYP");
    string LDA("LDA");
    transform(mtd.begin(), mtd.end(), mtd.begin(), ::toupper);
    if (mtd.compare(HF) == 0)
    {
        return 1;
    }
    else if (mtd.compare(B3LYP) == 0)
    {
        return 21;
    }
    else if (mtd.compare(LDA) == 0)
    {
        return 22;
    }
    else
    {
        cout << "Method not included or wrong" << endl;
        abort();
    }
}
