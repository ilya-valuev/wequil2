
#pragma hdrstop
#include <condefs.h>


//---------------------------------------------------------------------------
USEUNIT("interact.cpp");
USEUNIT("plasma.cpp");
USEUNIT("plasmar.cpp");
USEUNIT("..\myinc\record.cpp");
USEUNIT("..\myinc\common.cpp");
USEUNIT("props\anvdistr.cpp");
USEUNIT("props\anpair.cpp");
USEUNIT("props\ancurcur.cpp");
USEUNIT("props\andstruc.cpp");
USEUNIT("statist.cpp");
USEUNIT("..\myinc\GRAPH\Vector_3.cpp");
USEUNIT("..\myinc\RPARAM\RPARAM.C");
USEUNIT("..\myinc\erf.c");
USEUNIT("nanalyse.cpp");
USEUNIT("analyser.cpp");
USEUNIT("..\myinc\FOUR\four.c");
USEUNIT("props\antherm.cpp");
USEUNIT("props\anvcorr.cpp");
//---------------------------------------------------------------------------
extern int main_wan(int argc, char** argv);

#include <new.h>

#pragma argsused
int main(int argc, char **argv)
{
        set_new_handler(0);  // lets operator new return NULL in case of low memory  
        return main_wan(argc,argv);
}
