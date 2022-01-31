//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "ChildGrafica.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "ChildWin"
#pragma resource "*.dfm"
TMDIChild2 *MDIChild2;
//---------------------------------------------------------------------------
__fastcall TMDIChild2::TMDIChild2(TComponent* Owner)
        : TMDIChild(Owner)
{
}
//---------------------------------------------------------------------------
