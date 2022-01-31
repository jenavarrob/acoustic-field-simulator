//---------------------------------------------------------------------------
#ifndef ChildGraficaH
#define ChildGraficaH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "ChildWin.h"
#include <Chart.hpp>
#include <ExtCtrls.hpp>
#include <Series.hpp>
#include <TeEngine.hpp>
#include <TeeProcs.hpp>
//---------------------------------------------------------------------------
class TMDIChild2 : public TMDIChild
{
__published:	// IDE-managed Components
        TChart *Chart1;
        TLineSeries *Series1;
        TLineSeries *Series2;
        TLineSeries *Series3;
private:	// User declarations
public:		// User declarations
        __fastcall TMDIChild2(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TMDIChild2 *MDIChild2;
//---------------------------------------------------------------------------
#endif
