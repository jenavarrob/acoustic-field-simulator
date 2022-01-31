//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <dir.h>
#include <string.h>
#include "Ayuda.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TAyudaForm *AyudaForm;
//---------------------------------------------------------------------------
__fastcall TAyudaForm::TAyudaForm(TComponent* Owner)
        : TForm(Owner)
{
//HTML1->RequestDoc("file:d:\\Emeterio\\Pro-tec2\\Programas\\SiCampA\\Ayuda.html");
AnsiString Rotas;
char *pathe,*ruta;
char hoja[10];
char drive[MAXDRIVE];
char dir[MAXDIR];

pathe="file:";
strcpy(hoja,"Ayuda.html");
Rotas=Application->ExeName.c_str();
ruta=Rotas.c_str();
fnsplit(ruta, drive,dir, 0, 0);
StrCat(pathe,drive);
StrCat(pathe,dir);
StrCat(pathe,hoja);
HTML1->RequestDoc(pathe);
}
//---------------------------------------------------------------------------
