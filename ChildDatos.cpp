//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "ChildGrafica2D.h"
#include "ChildDatos.h"
#include "Main.h"
#include "ChildGrafica.h"
#include "Clases.h"
#include <math.h>
#include "ChildImagen.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "ChildWin"
#pragma resource "*.dfm"
TMDIChild1 *MDIChild1;
//---------------------------------------------------------------------------
__fastcall TMDIChild1::TMDIChild1(TComponent* Owner)
        : TMDIChild(Owner)
{
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::Button_LiquidosClick(TObject *Sender)
{
TMDIChild4 *ChildsContorno;
TMDIChild3 *ChildsPres2D;
TMDIChild2 *ChildsPresX,*ChildsPresY,*ChildsPresZ;
float **Real,**Imag,**FReal,**FImag,**GReal,**GImag,**FFR,**FFI,**Pres,*x,*y,*sy;
float A,delta_omega,delta_te,factor,r,z,ro,epsilon0,epsilon1,tei,a,b,densidad,paso_axial
,rango_axial_ini,rango_axial_fin,f,xrec,yrec,maximo=0,pos_maximo;
int i,j,k,rango,N,num_splines,pixels1,pixels2,pixels3,S;
FILE *APresion2D;

MainForm->StatusBar->SimpleText="Procesando";
//inicializa variables
r=StrToFloat(Edit_r->Text);
a=StrToFloat(Edit_a->Text);
b=StrToFloat(Edit_b->Text);
N=StrToInt(Edit_N->Text);
ro=StrToFloat(Edit_Ro->Text);
z=StrToFloat(Edit_z->Text);
rango_axial_ini=StrToFloat(Edit_zini->Text);
rango_axial_fin=StrToFloat(Edit_zfin->Text);
f=StrToFloat(Edit_f->Text);
xrec=StrToFloat(Edit_x->Text);
yrec=StrToFloat(Edit_y->Text);
Parametros(N,r,a,b,xrec,yrec,&A,&delta_te,&delta_omega,&rango,RG_Liquidos->ItemIndex);  //obten parametros necesarios
factor=A*15;
epsilon0=0.01;
epsilon1=0.01;
num_splines=10;
paso_axial=1;
//inicializa matrices y vectores
Real=matrix(-rango,rango-1,-rango,rango-1);
Imag=matrix(-rango,rango-1,-rango,rango-1);
FReal=matrix(-rango,rango-1,-rango,rango-1);
FImag=matrix(-rango,rango-1,-rango,rango-1);
GReal=matrix(-rango,rango-1,-rango,rango-1);
GImag=matrix(-rango,rango-1,-rango,rango-1);
FFR=matrix(-rango,rango-1,-rango,rango-1);
FFI=matrix(-rango,rango-1,-rango,rango-1);
Pres=matrix(-rango,rango-1,-rango,rango-1);
x=vector(1,num_splines);
y=vector(1,num_splines);
sy=vector(1,num_splines);  //vector donde estaran resultados de los splines

Genera_Splines(x,y,sy,num_splines,epsilon0,epsilon1);
Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Liquidos->ItemIndex);
Obten_TFuncion_GReen(GReal,GImag,delta_omega,rango,ro,z);
Copia_Ms(FReal,FImag,Real,Imag,rango);
TFourier(FReal,FImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,1);
Obten_Presion(Pres,FFR,FFI,FReal,FImag,GReal,GImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,rango);
//Inicializa Childs para graficar
if(CheckListBoxLiq->Checked[0]==true)  //presiones transversales
 {
    ChildsPresX=new TMDIChild2(Application);
    ChildsPresX->Caption="Presion Perfil Transversal X";
    ChildsPresY=new TMDIChild2(Application);
    ChildsPresY->Caption="Presion Perfil Transversal Y";
    for(i=-rango;i<=rango-1;i++)
     {
        tei=i*delta_te;
        ChildsPresX->Series1->AddXY(tei,Pres[i][0],"",clBlack);     //compresion
        ChildsPresY->Series1->AddXY(tei,Pres[0][i],"",clBlack);
     }
 }
if(CheckListBoxLiq->Checked[1]==true)// //Grafica presion acustica   3D
 {
    ChildsPres2D = new TMDIChild3(Application);
    ChildsPres2D->Caption="Presion Acustica";
    ChildsPres2D->Inicializa_Child(factor,rango,N);
    ChildsPres2D->Grafica3D(Pres,N,N);
 }
if(CheckListBoxLiq->Checked[2]==true)      //graficar perfil axial
 {
    ChildsPresZ=new TMDIChild2(Application);
    ChildsPresZ->Caption="Presion Perfil Axial";
    Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Liquidos->ItemIndex);
    for(k=rango_axial_ini;k<rango_axial_fin;k=k+paso_axial)
     {
       Copia_Ms(FReal,FImag,Real,Imag,rango);
       Obten_TFuncion_GReen(GReal,GImag,delta_omega,rango,ro,k);
       TFourier(FReal,FImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,1);
       Obten_Presion(Pres,FFR,FFI,FReal,FImag,GReal,GImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,rango);
       ChildsPresZ->Series1->AddXY(k,Pres[0][0],"",clBlack);
       if(Pres[0][0]>maximo)
        {
          maximo=Pres[0][0];
          pos_maximo=k;
        }
       Application->ProcessMessages();
     }
   for(i=ChildsPresZ->Series1->MinYValue();i<ChildsPresZ->Series1->MaxYValue();i++)
      ChildsPresZ->Series3->AddXY(pos_maximo,i,"",clRed);
 }
if(CheckListBoxLiq->Checked[3]==true)  //grafica presion contorno
{
  ChildsContorno=new TMDIChild4(Application);
  ChildsContorno->Caption="Presion Contorno";
  ChildsContorno->Width=rango_axial_fin-rango_axial_ini;
  ChildsContorno->Height=2*rango;
  Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Liquidos->ItemIndex);
  for(i=0,k=rango_axial_ini;k<rango_axial_fin;k=k+paso_axial,i++)
   {
     Copia_Ms(FReal,FImag,Real,Imag,rango);
     Obten_TFuncion_GReen(GReal,GImag,delta_omega,rango,ro,k);
     TFourier(FReal,FImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,1);
     Obten_Presion(Pres,FFR,FFI,FReal,FImag,GReal,GImag,N,r,a,b,xrec,yrec,RG_Liquidos->ItemIndex,rango);
     for(j=0;j<2*rango;j++)
           ChildsContorno->Image1->Canvas->Pixels[i][j]=TColor(RGB(Pres[j-rango][0]*255*(1/ro),Pres[j-rango][0]*255*(1/ro),Pres[j-rango][0]*255*(1/ro)));
     Application->ProcessMessages();
   }
}
if(CheckListBoxLiq->Checked[4]==true)  //guardar en archivo
 {
   if(MainForm->SaveDialog1->Execute())
    {
       APresion2D=fopen(MainForm->SaveDialog1->FileName.c_str(), "wt");
       Datos_Archivo(APresion2D,Pres,rango);
       fclose(APresion2D);
    }
 }
MainForm->StatusBar->SimpleText="Inactivo";

//libera memoria
free_matrix(Real,-rango,rango-1,-rango,rango-1);
free_matrix(Imag,-rango,rango-1,-rango,rango-1);
free_matrix(FReal,-rango,rango-1,-rango,rango-1);
free_matrix(FImag,-rango,rango-1,-rango,rango-1);
free_matrix(GReal,-rango,rango-1,-rango,rango-1);
free_matrix(GImag,-rango,rango-1,-rango,rango-1);
free_matrix(FFR,-rango,rango-1,-rango,rango-1);
free_matrix(FFI,-rango,rango-1,-rango,rango-1);
free_matrix(Pres,-rango,rango-1,-rango,rango-1);
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::Button2Click(TObject *Sender)
{
Close();
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::Button4Click(TObject *Sender)
{
Close();
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::Button_SolidosClick(TObject *Sender)
{
TMDIChild4 *ChildsContornoComp,*ChildsContornoCorte;
TMDIChild3 *ChildsPres2DCom,*ChildsPres2DCor;
TMDIChild2 *ChildsPresX,*ChildsPresY,*ChildsPresZ,*ChildsECoCo;
float **Real,**Imag,**FReal,**FImag,**FRealTemp,**FImagTemp,**GReal,**GImag,**GReal2,**GImag2,**FFR,**FFI,**Pres,**Pres1,*x,*y,*sy;
float A,delta_omega,delta_te,factor,r,z,ro,epsilon0,epsilon1,tei,a,b,alfa,beta,densidad,w_ini,w_fin,f
,itera_rango,paso_axial,rango_axial_ini,rango_axial_fin,paso_rango_w,xrec,yrec;
int i,j,k,rango,rango_axial,N,num_splines;
FILE *APresion2D;

MainForm->StatusBar->SimpleText="Procesando";
//inicializa variables
r=StrToFloat(Edit_rS->Text);
a=StrToFloat(Edit_aS->Text);
b=StrToFloat(Edit_bS->Text);
N=StrToInt(Edit_NS->Text);
z=StrToFloat(Edit_zS->Text);
f=StrToFloat(Edit_fS->Text);
rango_axial_ini=StrToFloat(Edit_ziniS->Text);
rango_axial_fin=StrToFloat(Edit_zfinS->Text);
alfa=StrToFloat(Edit_alfa->Text);
beta=StrToFloat(Edit_beta->Text);
densidad=StrToFloat(Edit_RoS->Text);
w_ini=StrToFloat(Edit_wini->Text);
w_fin=StrToFloat(Edit_wfin->Text);
xrec=StrToFloat(Edit_xS->Text);
yrec=StrToFloat(Edit_yS->Text);
Parametros(N,r,a,b,xrec,yrec,&A,&delta_te,&delta_omega,&rango,RG_Solidos->ItemIndex);  //obten parametros necesarios
factor=A*15;
epsilon0=0.01;
epsilon1=0.01;
num_splines=10;
paso_axial=1;
paso_rango_w=0.05;
//inicializa matrices y vectores
Real=matrix(-rango,rango-1,-rango,rango-1);
Imag=matrix(-rango,rango-1,-rango,rango-1);
FReal=matrix(-rango,rango-1,-rango,rango-1);
FImag=matrix(-rango,rango-1,-rango,rango-1);
FRealTemp=matrix(-rango,rango-1,-rango,rango-1);
FImagTemp=matrix(-rango,rango-1,-rango,rango-1);
GReal=matrix(-rango,rango-1,-rango,rango-1);
GImag=matrix(-rango,rango-1,-rango,rango-1);
GReal2=matrix(-rango,rango-1,-rango,rango-1);
GImag2=matrix(-rango,rango-1,-rango,rango-1);
FFR=matrix(-rango,rango-1,-rango,rango-1);
FFI=matrix(-rango,rango-1,-rango,rango-1);
Pres=matrix(-rango,rango-1,-rango,rango-1);
Pres1=matrix(-rango,rango-1,-rango,rango-1);
x=vector(1,num_splines);
y=vector(1,num_splines);
sy=vector(1,num_splines);  //vector donde estaran resultados de los splines

Genera_Splines(x,y,sy,num_splines,epsilon0,epsilon1);
Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Solidos->ItemIndex);
Obten_TFuncion_GReenSS(GReal,GImag,GReal2,GImag2,delta_omega,rango,densidad,z,alfa,beta,f,0);
Copia_Ms(FRealTemp,FImagTemp,Real,Imag,rango);
TFourier(FRealTemp,FImagTemp,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,1);
Obten_Presion(Pres,FFR,FFI,FRealTemp,FImagTemp,GReal,GImag,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);  //compresion
Obten_Presion(Pres1,FFR,FFI,FRealTemp,FImagTemp,GReal2,GImag2,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango); //corte
//Inicializa Childs para graficar
if(CheckListBoxSol->Checked[0]==true)      //grafica presion acustica secciones transversales
 {
    ChildsPresX=new TMDIChild2(Application);
    ChildsPresX->Caption="Desplazamiento Perfil Transversal X";
    ChildsPresY=new TMDIChild2(Application);
    ChildsPresY->Caption="Desplazamiento Perfil Transversal Y";
    for(i=-rango;i<=rango-1;i++)
     {
       tei=i*delta_te;
       ChildsPresX->Series2->AddXY(tei,Pres1[i][0],"",clBlue);  //corte
       ChildsPresY->Series2->AddXY(tei,Pres1[0][i],"",clBlue);
       ChildsPresX->Series1->AddXY(tei,Pres[i][0],"",clBlack);  //compresion
       ChildsPresY->Series1->AddXY(tei,Pres[0][i],"",clBlack);
     }
 }
if(CheckListBoxSol->Checked[1]==true)      //Grafica presion acustica   3D
 {
    ChildsPres2DCom = new TMDIChild3(Application);
    ChildsPres2DCom->Caption="Desplazamiento Acustico Compresion";
    ChildsPres2DCom->Inicializa_Child(factor,rango,N);
    ChildsPres2DCom->Grafica3D(Pres,N,N);
    ChildsPres2DCor = new TMDIChild3(Application);
    ChildsPres2DCor->Caption="Desplazamiento Acustico Corte";
    ChildsPres2DCor->Inicializa_Child(factor,rango,N);
    ChildsPres2DCor->Grafica3D(Pres1,N,N);
 }
if(CheckListBoxSol->Checked[2]==true)      //graficar perfil axial
 {
    ChildsPresZ=new TMDIChild2(Application);
    ChildsPresZ->Caption="Desplazamiento Perfil Axial";
    Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Solidos->ItemIndex);
    for(k=rango_axial_ini;k<rango_axial_fin;k=k+paso_axial)
     {
        Copia_Ms(FReal,FImag,Real,Imag,rango);
        Obten_TFuncion_GReenSS(GReal,GImag,GReal2,GImag2,delta_omega,rango,densidad,k,alfa,beta,f,0);
        Copia_Ms(FRealTemp,FImagTemp,Real,Imag,rango);
        TFourier(FRealTemp,FImagTemp,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,1);
        Obten_Presion(Pres,FFR,FFI,FRealTemp,FImagTemp,GReal,GImag,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
        Obten_Presion(Pres1,FFR,FFI,FRealTemp,FImagTemp,GReal2,GImag2,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
        ChildsPresZ->Series1->AddXY(k,Pres[0][0],"",clBlack);    //compresion
        ChildsPresZ->Series2->AddXY(k,Pres1[0][0],"",clBlue);   //corte
        Application->ProcessMessages();
     }
 }
if(CheckListBoxSol->Checked[3]==true)   //graficar contorno
 {
  ChildsContornoComp=new TMDIChild4(Application);
  ChildsContornoComp->Caption="Desplazamiento-Compresion";
  ChildsContornoCorte=new TMDIChild4(Application);
  ChildsContornoCorte->Caption="Desplazamiento-Corte";
  ChildsContornoComp->Width=rango_axial_fin-rango_axial_ini;
  ChildsContornoComp->Height=2*rango;
  ChildsContornoCorte->Width=rango_axial_fin-rango_axial_ini;
  ChildsContornoCorte->Height=2*rango;
  Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Liquidos->ItemIndex);
  for(i=0,k=rango_axial_ini;k<rango_axial_fin;k=k+paso_axial,i++)
   {
     Copia_Ms(FReal,FImag,Real,Imag,rango);
     Obten_TFuncion_GReenSS(GReal,GImag,GReal2,GImag2,delta_omega,rango,densidad,k,alfa,beta,f,0);
     Copia_Ms(FRealTemp,FImagTemp,Real,Imag,rango);
     TFourier(FRealTemp,FImagTemp,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,1);
     Obten_Presion(Pres,FFR,FFI,FRealTemp,FImagTemp,GReal,GImag,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
     Obten_Presion(Pres1,FFR,FFI,FRealTemp,FImagTemp,GReal2,GImag2,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
     for(j=0;j<2*rango;j++)
        {
           ChildsContornoComp->Image1->Canvas->Pixels[i][j]=TColor(RGB(Pres[j-rango][0]*1000*(1/densidad),Pres[j-rango][0]*1000*(1/densidad),Pres[j-rango][0]*1000*(1/densidad)));
           ChildsContornoCorte->Image1->Canvas->Pixels[i][j]=TColor(RGB(Pres1[j-rango][0]*1000*(1/densidad),Pres1[j-rango][0]*1000*(1/densidad),Pres1[j-rango][0]*1000*(1/densidad)));
        }
     Application->ProcessMessages();
   }
 }
if(CheckListBoxSol->Checked[4]==true)  //guardar en archivo
 {
   if(MainForm->SaveDialog1->Execute())
    {
       APresion2D=fopen(MainForm->SaveDialog1->FileName.c_str(), "wt");
       Datos_Archivo(APresion2D,Pres,rango);
       fclose(APresion2D);
    }
 }
if(CheckListBoxSol->Checked[5]==true)  //Espectro de Onda Compresion-Corte
{
float tei,tej,suma1,suma2;
  ChildsECoCo=new TMDIChild2(Application);
  ChildsECoCo->Caption="Espectro de Onda Compresion-Corte";
  Obten_Funcion_V(Real,Imag,delta_te,epsilon0,epsilon1,num_splines,r,rango,sy,x,y,a,b,xrec,yrec,RG_Solidos->ItemIndex);
    for(itera_rango=w_ini;itera_rango<w_fin;itera_rango=itera_rango+paso_rango_w)
     {
        Copia_Ms(FReal,FImag,Real,Imag,rango);
        Copia_Ms(FRealTemp,FImagTemp,Real,Imag,rango);
        Obten_TFuncion_GReenSS(GReal,GImag,GReal2,GImag2,delta_omega,rango,densidad,z,alfa,beta,itera_rango,1);
        TFourier(FReal,FImag,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,1);
        Obten_Presion(Pres,FFR,FFI,FReal,FImag,GReal,GImag,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
        TFourier(FRealTemp,FImagTemp,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,1);
        Obten_Presion(Pres1,FFR,FFI,FRealTemp,FImagTemp,GReal2,GImag2,N,r,a,b,xrec,yrec,RG_Solidos->ItemIndex,rango);
        //Para obtener Espectro, integral sobre el resultado de la presion
        suma1=0;
        suma2=0;
        for(i=-rango;i<=rango-1;i++)
          for(j=-rango;j<=rango-1;j++)
            {
               if(RG_Solidos->ItemIndex==0)  //circulares
                 {
                  tei=i*delta_te/r;
                  tej=j*delta_te/r;
                  if((tei*tei)+(tej*tej)<1)
                    {
                    suma1+=(Pres[i][j]*delta_te*delta_te)/(r*r);
                    suma2+=(Pres1[i][j]*delta_te*delta_te)/(r*r);
                    }
                 }
               if(RG_Solidos->ItemIndex==1)  //elipticos
                 {
                  tei=i*delta_te/a;
                  tej=j*delta_te/b;
                  if((tei*tei)+(tej*tej)<1)
                   {
                    suma1+=(Pres[i][j]*delta_te*delta_te)/(a*b);
                    suma2+=(Pres1[i][j]*delta_te*delta_te)/(a*b);
                   }
                 }
            }
         suma1=suma1/M_PI;
         suma2=suma2/M_PI;
         ChildsECoCo->Series1->AddXY(itera_rango,suma1,"",clBlack);  //compresion
         ChildsECoCo->Series2->AddXY(itera_rango,suma2,"",clBlue);   //corte
         Application->ProcessMessages();
    }
}

MainForm->StatusBar->SimpleText="Inactivo";

//libera memoria
free_matrix(Real,-rango,rango-1,-rango,rango-1);
free_matrix(Imag,-rango,rango-1,-rango,rango-1);
free_matrix(FReal,-rango,rango-1,-rango,rango-1);
free_matrix(FImag,-rango,rango-1,-rango,rango-1);
free_matrix(GReal,-rango,rango-1,-rango,rango-1);
free_matrix(GImag,-rango,rango-1,-rango,rango-1);
free_matrix(GReal2,-rango,rango-1,-rango,rango-1);
free_matrix(GImag2,-rango,rango-1,-rango,rango-1);
free_matrix(FFR,-rango,rango-1,-rango,rango-1);
free_matrix(FFI,-rango,rango-1,-rango,rango-1);
free_matrix(Pres,-rango,rango-1,-rango,rango-1);
free_matrix(Pres1,-rango,rango-1,-rango,rango-1);
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::RG_LiquidosClick(TObject *Sender)
{
switch(RG_Liquidos->ItemIndex)
 {
   case 0:   //circular
         Edit_r->Enabled=true;
         Edit_a->Enabled=false;
         Edit_b->Enabled=false;
         Edit_x->Enabled=false;
         Edit_y->Enabled=false;
        break;
   case 1:   //eliptico
         Edit_r->Enabled=false;
         Edit_a->Enabled=true;
         Edit_b->Enabled=true;
         Edit_x->Enabled=false;
         Edit_y->Enabled=false;
        break;
   case 2:   //rectangular
         Edit_r->Enabled=false;
         Edit_a->Enabled=false;
         Edit_b->Enabled=false;
         Edit_x->Enabled=true;
         Edit_y->Enabled=true;
        break;
   default:
        break;
 }
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::CheckListBoxLiqClick(TObject *Sender)
{
int i;
bool habilitado;
for(i=0;i<5;i++)
 {
 habilitado=CheckListBoxLiq->Checked[i];
  switch(i)
   {
      case 0:        //Presion transversal 2D
             if(CheckListBoxLiq->Checked[1]==true || CheckListBoxLiq->Checked[4]==true)
                  break;
             else
                  Edit_z->Enabled=habilitado;
             break;
      case 1:        //Presion transversal 3D
             if(CheckListBoxLiq->Checked[0]==true || CheckListBoxLiq->Checked[4]==true)
                  break;
             else
                  Edit_z->Enabled=habilitado;
             break;
      case 2:        //perfil axial
             if(CheckListBoxLiq->Checked[3]==true)
                break;
             else
               {
                 Edit_zini->Enabled=habilitado;
                 Edit_zfin->Enabled=habilitado;
               }
             break;
      case 3:        //contorno
             if(CheckListBoxLiq->Checked[2]==true)
                break;
             else
               {
                 Edit_zini->Enabled=habilitado;
                 Edit_zfin->Enabled=habilitado;
               }
             break;
      case 4:        //guardar datos
             if(CheckListBoxLiq->Checked[0]==true || CheckListBoxLiq->Checked[1]==true)
                  break;
             else
                  Edit_z->Enabled=habilitado;
             break;
      default:
             break;

   }
 }
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::RG_SolidosClick(TObject *Sender)
{
switch(RG_Solidos->ItemIndex)
 {
   case 0:   //circular
         Edit_rS->Enabled=true;
         Edit_aS->Enabled=false;
         Edit_bS->Enabled=false;
         Edit_xS->Enabled=false;
         Edit_yS->Enabled=false;
        break;
   case 1:   //eliptico
         Edit_rS->Enabled=false;
         Edit_aS->Enabled=true;
         Edit_bS->Enabled=true;
         Edit_xS->Enabled=false;
         Edit_yS->Enabled=false;
        break;
   case 2:   //rectangular
         Edit_rS->Enabled=false;
         Edit_aS->Enabled=false;
         Edit_bS->Enabled=false;
         Edit_xS->Enabled=true;
         Edit_yS->Enabled=true;
        break;
   default:
        break;
 }
}
//---------------------------------------------------------------------------

void __fastcall TMDIChild1::CheckListBoxSolClick(TObject *Sender)
{
int i;
bool habilitado;
for(i=0;i<6;i++)
 {
 habilitado=CheckListBoxSol->Checked[i];
  switch(i)
   {
      case 0:        //Presion transversal 2D
             if(CheckListBoxSol->Checked[1]==true || CheckListBoxSol->Checked[4]==true)
                  break;
             else
                  Edit_zS->Enabled=habilitado;
             break;
      case 1:        //Presion transversal 3D
             if(CheckListBoxSol->Checked[0]==true || CheckListBoxSol->Checked[4]==true)
                  break;
             else
                  Edit_zS->Enabled=habilitado;
             break;
      case 2:        //perfil axial
             if(CheckListBoxSol->Checked[3]==true)
                break;
             else
               {
                 Edit_ziniS->Enabled=habilitado;
                 Edit_zfinS->Enabled=habilitado;
               }
             break;
      case 3:        //contorno
             if(CheckListBoxSol->Checked[2]==true)
                break;
             else
               {
                 Edit_ziniS->Enabled=habilitado;
                 Edit_zfinS->Enabled=habilitado;
               }
             break;
      case 4:        //guardar datos
             if(CheckListBoxSol->Checked[0]==true || CheckListBoxSol->Checked[1]==true)
                  break;
             else
                  Edit_zS->Enabled=habilitado;
             break;
      case 5:        //espectro de onda
             Edit_wini->Enabled=habilitado;
             Edit_wfin->Enabled=habilitado;
             break;
      default:
             break;

   }
 }
}
//---------------------------------------------------------------------------

