//---------------------------------------------------------------------------
#ifndef ChildGrafica2DH
#define ChildGrafica2DH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "ChildWin.h"
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TMDIChild3 : public TMDIChild
{
__published:	// IDE-managed Components
        TImage *Image1;
        void __fastcall Image1MouseDown(TObject *Sender,
          TMouseButton Button, TShiftState Shift, int X, int Y);
private:	// User declarations
public:		// User declarations
        __fastcall TMDIChild3(TComponent* Owner);
        void Grafica3D(float **F,int nr,int nc);
        void ReGrafica(void);
        float **Superficie;
        float fact;
        int N;
        void __fastcall Inicializa_Child(int factor,int rango,int Ene);
};
//---------------------------------------------------------------------------
extern PACKAGE TMDIChild3 *MDIChild3;
//---------------------------------------------------------------------------
#endif
