//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "Main.h"
#include "About.h"
#include "Clases.h"
#include "ChildDatos.h"
#include "ChildGrafica.h"
#include "Ayuda.h"
//---------------------------------------------------------------------------
#pragma resource "*.dfm"
TMainForm *MainForm;
//---------------------------------------------------------------------------

__fastcall TMainForm::TMainForm(TComponent *Owner)
	: TForm(Owner)
{
        CopyItem->Enabled=true;
        ToolButton5->Enabled=true;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CreateMDIChild(String Name)
{
TMDIChild *Child1;
Child1 = new TMDIChild1(Application);
Child1->Caption = Name;
Child1->AutoSize=true;
}
//---------------------------------------------------------------------------
void __fastcall TMainForm::CreateMDIChild1(String Name)
{
FILE *ADatos;
TMDIChild1 *Child1;
AnsiString linea;
char val;
  //--- create a new MDI child window ----
  Child1 = new TMDIChild1(Application);
  Child1->Caption = Name;
  Child1->AutoSize=true;
  if (FileExists (Name))
    {
        ADatos=fopen(Name.c_str(), "rt");
        while(!feof(ADatos))
          {
            fscanf(ADatos, "%c",&val);
            linea="";
            switch(val)
              {
                case 'N':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_N->Text=linea;
                      break;
                case 'z':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_z->Text=linea;
                      break;
                case 'i':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_zini->Text=linea;
                      break;
                case 'f':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_zfin->Text=linea;
                      break;
                case 'R':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      Child1->RG_Liquidos->ItemIndex=StrToInt(val);
                      break;
                case 'D':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_Ro->Text=linea;
                      break;
                case 'x':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_f->Text=linea;
                      break;
                case 'C':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_r->Text=linea;
                      break;
                case 'A':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_a->Text=linea;
                      break;
                case 'B':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_b->Text=linea;
                      break;
                case '#':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_x->Text=linea;
                      break;
                case '$':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_y->Text=linea;
                      break;

                case '0':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxLiq->Checked[0]=true;
                      else
                        Child1->CheckListBoxLiq->Checked[0]=false;
                      break;
                case '1':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxLiq->Checked[1]=true;
                      else
                        Child1->CheckListBoxLiq->Checked[1]=false;
                      break;
                case '2':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxLiq->Checked[2]=true;
                      else
                        Child1->CheckListBoxLiq->Checked[2]=false;
                      break;
                case '3':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxLiq->Checked[3]=true;
                      else
                        Child1->CheckListBoxLiq->Checked[3]=false;
                      break;
                case '4':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxLiq->Checked[4]=true;
                      else
                        Child1->CheckListBoxLiq->Checked[4]=false;
                      break;
//solidos
                case 'M':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_NS->Text=linea;
                      break;
                case 'e':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_zS->Text=linea;
                      break;
                case 'g':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_ziniS->Text=linea;
                      break;
                case 'h':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_zfinS->Text=linea;
                      break;
                case 'j':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      Child1->RG_Solidos->ItemIndex=StrToInt(val);
                      break;
                case 'k':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_RoS->Text=linea;
                      break;
                case 'l':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_fS->Text=linea;
                      break;
                case 'o':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_rS->Text=linea;
                      break;
                case 'p':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_aS->Text=linea;
                      break;
                case 'q':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_bS->Text=linea;
                      break;
                case '%':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_xS->Text=linea;
                      break;
                case '&':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_yS->Text=linea;
                      break;

                case '5':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[0]=true;
                      else
                        Child1->CheckListBoxSol->Checked[0]=false;
                      break;
                case '6':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[1]=true;
                      else
                        Child1->CheckListBoxSol->Checked[1]=false;
                      break;
                case '7':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[2]=true;
                      else
                        Child1->CheckListBoxSol->Checked[2]=false;
                      break;
                case '8':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[3]=true;
                      else
                        Child1->CheckListBoxSol->Checked[3]=false;
                      break;
                case '9':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[4]=true;
                      else
                        Child1->CheckListBoxSol->Checked[4]=false;
                      break;
                case 's':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      if(StrToInt(val)==1)
                        Child1->CheckListBoxSol->Checked[5]=true;
                      else
                        Child1->CheckListBoxSol->Checked[5]=false;
                      break;
                case 'u':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_alfa->Text=linea;
                      break;
                case 'v':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_beta->Text=linea;
                      break;
                case 'w':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_wini->Text=linea;
                      break;
                case 'Y':
                      fscanf(ADatos, "%c",&val);
                      fscanf(ADatos, "%c",&val);
                      while(val!='\n')
                       {
                         linea+=val;
                         fscanf(ADatos, "%c",&val);
                       }
                      Child1->Edit_wfin->Text=linea;
                      break;
                default: break;
              }
          }
       fclose(ADatos);
    }
  else
    ShowMessage("Archivo no existe");
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FileNew1Execute(TObject *Sender)
{
	CreateMDIChild("Datos" + IntToStr(MDIChildCount + 1));
        FileSaveItem->Enabled=true;
        ToolButton2->Enabled=true;
        CopyItem->Enabled=true;
        ToolButton5->Enabled=true;

}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FileOpen1Execute(TObject *Sender)
{
	if (OpenDialog->Execute())
		CreateMDIChild1(OpenDialog->FileName);
        FileSaveItem->Enabled=true;
        ToolButton2->Enabled=true;
        CopyItem->Enabled=true;
        ToolButton5->Enabled=true;
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::HelpAbout1Execute(TObject *Sender)
{
	AboutBox->ShowModal();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FileExit1Execute(TObject *Sender)
{
	Close();
}
//---------------------------------------------------------------------------

void __fastcall CreateMDIChild2(String Name)
{
        TMDIChild2 *Child2;
	//--- create a new MDI child window ----
	Child2 = new TMDIChild2(Application);
	Child2->Caption = Name;
        Child2->AutoSize=true;
}

void __fastcall TMainForm::CerrarTodo1Click(TObject *Sender)
{
 for(int i = 0; i < MDIChildCount; i++)
    MDIChildren[i]->Close();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::FileSaveItemClick(TObject *Sender)
{
FILE *AConfig;
TMDIChild1 *Child1;
TMDIChild *Child;
AnsiString linea,linea1;
Graphics::TBitmap *Bitmap1 = new Graphics::TBitmap();
float val;
Child1=(TMDIChild1 *)ActiveMDIChild;
linea="TMDIChild1";

if(linea==String(Child1->ClassName()))
  {
    if (SaveDialog1->Execute())
     {
        AConfig=fopen(SaveDialog1->FileName.c_str(), "wt");
        fprintf(AConfig,"N=");
        val=StrToFloat(Child1->Edit_N->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"z=");
        val=StrToFloat(Child1->Edit_z->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"i=");
        val=StrToFloat(Child1->Edit_zini->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"f=");
        val=StrToFloat(Child1->Edit_zfin->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"R=");
        val=StrToFloat(Child1->RG_Liquidos->ItemIndex);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"D=");
        val=StrToFloat(Child1->Edit_Ro->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"x=");
        val=StrToFloat(Child1->Edit_f->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"C=");
        val=StrToFloat(Child1->Edit_r->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"A=");
        val=StrToFloat(Child1->Edit_a->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"B=");
        val=StrToFloat(Child1->Edit_b->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"#=");
        val=StrToFloat(Child1->Edit_x->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"$=");
        val=StrToFloat(Child1->Edit_y->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"0=");
        val=Child1->CheckListBoxLiq->Checked[0];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"1=");
        val=Child1->CheckListBoxLiq->Checked[1];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"2=");
        val=Child1->CheckListBoxLiq->Checked[2];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"3=");
        val=Child1->CheckListBoxLiq->Checked[3];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"4=");
        val=Child1->CheckListBoxLiq->Checked[4];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"M=");
        val=StrToFloat(Child1->Edit_NS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"e=");
        val=StrToFloat(Child1->Edit_zS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"g=");
        val=StrToFloat(Child1->Edit_ziniS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"h=");
        val=StrToFloat(Child1->Edit_zfinS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"j=");
        val=StrToFloat(Child1->RG_Solidos->ItemIndex);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"k=");
        val=StrToFloat(Child1->Edit_RoS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"l=");
        val=StrToFloat(Child1->Edit_fS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"o=");
        val=StrToFloat(Child1->Edit_rS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"p=");
        val=StrToFloat(Child1->Edit_aS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"q=");
        val=StrToFloat(Child1->Edit_bS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"%=");
        val=StrToFloat(Child1->Edit_xS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"&=");
        val=StrToFloat(Child1->Edit_yS->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"5=");
        val=Child1->CheckListBoxSol->Checked[0];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"6=");
        val=Child1->CheckListBoxSol->Checked[1];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"7=");
        val=Child1->CheckListBoxSol->Checked[2];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"8=");
        val=Child1->CheckListBoxSol->Checked[3];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"9=");
        val=Child1->CheckListBoxSol->Checked[4];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"s=");
        val=Child1->CheckListBoxSol->Checked[5];
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"u=");
        val=StrToFloat(Child1->Edit_alfa->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"v=");
        val=StrToFloat(Child1->Edit_beta->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"w=");
        val=StrToFloat(Child1->Edit_wini->Text);
        fprintf(AConfig,"%f\n",val);
        fprintf(AConfig,"Y=");
        val=StrToFloat(Child1->Edit_wfin->Text);
        fprintf(AConfig,"%f\n",val);
       fclose(AConfig);
     }
  }
Child=(TMDIChild *)ActiveMDIChild;
linea="TMDIChild2";
linea1="TMDIChild3";
if(linea==String(Child->ClassName()) ||linea1==String(Child->ClassName()))
  {
    if (SavePictureDialog1->Execute())
      {
         Bitmap1=Child->GetFormImage();
         Bitmap1->SaveToFile(SavePictureDialog1->FileName);
      }
  }
}
//---------------------------------------------------------------------------



void __fastcall TMainForm::Contenido1Click(TObject *Sender)
{
  AyudaForm->ShowModal();
}
//---------------------------------------------------------------------------

void __fastcall TMainForm::CopyItemClick(TObject *Sender)
{
  Clipboard()->Assign(ActiveMDIChild->GetFormImage());
}
//---------------------------------------------------------------------------


