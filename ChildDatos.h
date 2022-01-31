//---------------------------------------------------------------------------
#ifndef ChildDatosH
#define ChildDatosH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include "ChildWin.h"
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <CheckLst.hpp>
#include <Tabnotbk.hpp>
//---------------------------------------------------------------------------
class TMDIChild1 : public TMDIChild
{
__published:	// IDE-managed Components
        TTabbedNotebook *TNotebook;
        TLabel *Label22;
        TLabel *Label23;
        TLabel *Labeld;
        TLabel *Label24;
        TLabel *Label25;
        TLabel *Label26;
        TLabel *Label27;
        TLabel *Label28;
        TLabel *Label29;
        TRadioGroup *RG_Liquidos;
        TBevel *Bevel3;
        TBevel *Bevel4;
        TEdit *Edit_N;
        TEdit *Edit_z;
        TEdit *Edit_zini;
        TEdit *Edit_Ro;
        TEdit *Edit_zfin;
        TEdit *Edit_f;
        TLabel *Label30;
        TEdit *Edit_r;
        TLabel *Label31;
        TEdit *Edit_a;
        TLabel *Label32;
        TEdit *Edit_b;
        TLabel *Label33;
        TEdit *Edit_x;
        TLabel *Label34;
        TEdit *Edit_y;
        TBevel *Bevel5;
        TCheckListBox *CheckListBoxLiq;
        TLabel *Label35;
        TButton *Button_Liquidos;
        TButton *Button2;
        TLabel *Label36;
        TLabel *Label37;
        TLabel *Label38;
        TLabel *Label39;
        TLabel *Label40;
        TLabel *Label41;
        TLabel *Label42;
        TLabel *Label43;
        TLabel *Label44;
        TBevel *Bevel6;
        TBevel *Bevel7;
        TLabel *Label45;
        TLabel *Label46;
        TLabel *Label47;
        TLabel *Label48;
        TLabel *Label49;
        TRadioGroup *RG_Solidos;
        TEdit *Edit_NS;
        TEdit *Edit_zS;
        TEdit *Edit_ziniS;
        TEdit *Edit_RoS;
        TEdit *Edit_zfinS;
        TEdit *Edit_fS;
        TEdit *Edit_rS;
        TEdit *Edit_aS;
        TEdit *Edit_bS;
        TEdit *Edit_xS;
        TEdit *Edit_yS;
        TBevel *Bevel8;
        TCheckListBox *CheckListBoxSol;
        TLabel *Label50;
        TButton *Button_Solidos;
        TButton *Button4;
        TLabel *Label51;
        TLabel *Label52;
        TLabel *Label53;
        TEdit *Edit_wini;
        TEdit *Edit_wfin;
        TLabel *Label54;
        TLabel *Label55;
        TEdit *Edit_alfa;
        TEdit *Edit_beta;
        TBevel *Bevel9;
        TBevel *Bevel10;
        void __fastcall Button_LiquidosClick(TObject *Sender);
        void __fastcall Button2Click(TObject *Sender);
        void __fastcall Button4Click(TObject *Sender);
        void __fastcall Button_SolidosClick(TObject *Sender);
        void __fastcall RG_LiquidosClick(TObject *Sender);
        void __fastcall CheckListBoxLiqClick(TObject *Sender);
        void __fastcall RG_SolidosClick(TObject *Sender);
        void __fastcall CheckListBoxSolClick(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TMDIChild1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TMDIChild1 *MDIChild1;
//---------------------------------------------------------------------------
#endif
