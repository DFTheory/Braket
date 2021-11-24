#include "DFT.hpp"
using namespace std;

void RunMission(int mission, int method, int basis, double E, vector<Atom> inputGroup, int Electrons, string FN)
{
    if (mission == 1)
    {
        if (method == 1)
        {
            if (basis == 1)
            {
                //cout << "Lets see if it can get here :v" << endl;
                //cout << "This is FN: " << FN << endl;
                HartreeFock(E, inputGroup, Electrons, FN);
            }
        }
        else if (method > 20)
        {
            if (basis == 1)
            {
                DFT(E, inputGroup, Electrons, FN, method);
            }
        }
    }
}