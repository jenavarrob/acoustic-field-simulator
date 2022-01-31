//---------------------------------------------------------------------------
#ifndef AyudaH
#define AyudaH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <NMHTML.hpp>
#include <OleCtrls.hpp>
//---------------------------------------------------------------------------
class TAyudaForm : public TForm
{
__published:	// IDE-managed Components
        THTML *HTML1;
private:	// User declarations
public:		// User declarations
        __fastcall TAyudaForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TAyudaForm *AyudaForm;
//---------------------------------------------------------------------------
#endif
