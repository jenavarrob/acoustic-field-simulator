//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "ChildImagen.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "ChildWin"
#pragma resource "*.dfm"
TMDIChild4 *MDIChild4;
//---------------------------------------------------------------------------
__fastcall TMDIChild4::TMDIChild4(TComponent* Owner)
        : TMDIChild(Owner)
{
Bitmap = new Graphics::TBitmap;
}
//---------------------------------------------------------------------------
