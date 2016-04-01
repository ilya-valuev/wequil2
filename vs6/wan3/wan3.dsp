# Microsoft Developer Studio Project File - Name="wan3" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=wan3 - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "wan3.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "wan3.mak" CFG="wan3 - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "wan3 - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "wan3 - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "wan3 - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "../" /I "../../../ivutils/include" /I "../../../gridmd/include" /I "../../../systems/interact" /I "../../../extlib/include" /I "../../../fitting/include" /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "NO_BFSEEK" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "wan3 - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
F90=df.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "../" /I "../../../ivutils/include" /I "../../../gridmd/include" /I "../../../systems/interact" /I "../../../extlib/include" /I "../../../fitting/include" /D "_DEBUG" /D "WIN32" /D "_CONSOLE" /D "_MBCS" /D "NO_BFSEEK" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "wan3 - Win32 Release"
# Name "wan3 - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\analyser.cpp
# End Source File
# Begin Source File

SOURCE=..\..\analyser.h
# End Source File
# Begin Source File

SOURCE=..\..\props\ancurcur.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\ancurcur.h
# End Source File
# Begin Source File

SOURCE=..\..\props\andstruc.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\andstruc.h
# End Source File
# Begin Source File

SOURCE=..\..\props\anmfield.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\anmfield.h
# End Source File
# Begin Source File

SOURCE=..\..\props\anpair.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\anpair.h
# End Source File
# Begin Source File

SOURCE=..\..\props\antherm.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\antherm.h
# End Source File
# Begin Source File

SOURCE=..\..\props\antrjlog.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\antrjlog.h
# End Source File
# Begin Source File

SOURCE=..\..\props\anvcorr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\anvcorr.h
# End Source File
# Begin Source File

SOURCE=..\..\props\anvdistr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\props\anvdistr.h
# End Source File
# Begin Source File

SOURCE=..\..\interact.cpp
# End Source File
# Begin Source File

SOURCE=..\..\interact.h
# End Source File
# Begin Source File

SOURCE=..\..\nanalyse.cpp
# End Source File
# Begin Source File

SOURCE=..\..\pcycle.cpp
# End Source File
# Begin Source File

SOURCE=..\..\pcycle.h
# End Source File
# Begin Source File

SOURCE=..\..\plasma.cpp
# End Source File
# Begin Source File

SOURCE=..\..\plasma.h
# End Source File
# Begin Source File

SOURCE=..\..\plasmar.cpp
# End Source File
# Begin Source File

SOURCE=..\..\plasmar.h
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Group "ivutils"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\..\ivutils\src\common.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\common.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\erf_namd.c
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\four\four.c
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\four.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\listiter.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\pairhash.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\record.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\record.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\rparam\rparam.c
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\rparam.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\statist.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\statist.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\src\graph\vector_3.cpp
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\vector_3.h
# End Source File
# Begin Source File

SOURCE=..\..\..\ivutils\include\vector_n.h
# End Source File
# End Group
# Begin Source File

SOURCE=..\..\examples\analyse_pair.cfg
# End Source File
# End Target
# End Project
