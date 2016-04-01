
extern int main_weq(int argc, char** argv);
#ifdef BCC
#pragma hdrstop
#include <condefs.h>


//---------------------------------------------------------------------------
USEUNIT("interact.cpp");
USEUNIT("pcycle.cpp");
USEUNIT("plasma.cpp");
USEUNIT("plasmar.cpp");
USEUNIT("cequilr.cpp");
USEUNIT("..\myinc\record.cpp");
USEUNIT("..\myinc\common.cpp");
USEUNIT("..\myinc\GRAPH\Vector_3.cpp");
USEUNIT("..\myinc\RPARAM\RPARAM.C");
USEUNIT("analyse.cpp");
USEUNIT("..\myinc\erf.c");
USEUNIT("statist.cpp");
//---------------------------------------------------------------------------


#pragma argsused
# endif
int main(int argc, char **argv)
{
        return main_weq(argc,argv);
}


