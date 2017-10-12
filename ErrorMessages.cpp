
#include "ErrorMessages.H"

void LongLine()
{
    cout << "-----------------------------------------------------" << endl;
};

void Waterloo(string f, string err)
{
    LongLine();
    cout << "                    WATERLOO                         " << endl;
    LongLine();
    cout << "In function: " << f << endl;
    cout << "Message: " << err << endl;
};
