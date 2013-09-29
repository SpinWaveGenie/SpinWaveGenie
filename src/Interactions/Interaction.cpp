#include "Interaction.h"

using namespace std;


bool Interaction::operator<(const Interaction& other) const
{
    vector<string> sl1 = this->sublattices();
    vector<string> sl2 = other.sublattices();
    if (sl1.size() <= sl2.size())
    {
        for (size_t index = 0;index < sl1.size();index++)
        {
            if (sl1[index] < sl2[index])
                return true;
            else if (sl1[index] > sl2[index])
                return false;
        }
        return true;
    }
    else
    {
        for (size_t index = 0;index < sl2.size();index++)
        {
            if (sl1[index] < sl2[index])
                return true;
            else if (sl1[index] > sl2[index])
                return false;
        }
        return false;
    }
}