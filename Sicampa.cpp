//---------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
//---------------------------------------------------------------------
USEFORM("Main.cpp", MainForm);
USEFORM("ChildWin.cpp", MDIChild);
USERES("SiCampA.res");
USEFORM("about.cpp", AboutBox);
USEUNIT("Clases.cpp");
USEFORM("ChildDatos.cpp", MDIChild1);
USEFORM("ChildGrafica.cpp", MDIChild2);
USEFORM("ChildGrafica2D.cpp", MDIChild3);
USEFORM("ChildImagen.cpp", MDIChild4);
USEFORM("Ayuda.cpp", AyudaForm);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
	Application->Initialize();
	Application->CreateForm(__classid(TMainForm), &MainForm);
                 Application->CreateForm(__classid(TAboutBox), &AboutBox);
                 Application->CreateForm(__classid(TAyudaForm), &AyudaForm);
                 Application->Run();

	return 0;
}
//---------------------------------------------------------------------
