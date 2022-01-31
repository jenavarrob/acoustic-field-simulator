//---------------------------------------------------------------------------
#ifndef ChildImagenH
#define ChildImagenH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "ChildWin.h"
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TMDIChild4 : public TMDIChild
{
__published:	// IDE-managed Components
        TImage *Image1;
private:	// User declarations
public:		// User declarations
        __fastcall TMDIChild4(TComponent* Owner);
        Graphics::TBitmap *Bitmap;        
};
//---------------------------------------------------------------------------
extern PACKAGE TMDIChild4 *MDIChild4;
//---------------------------------------------------------------------------
#endif
