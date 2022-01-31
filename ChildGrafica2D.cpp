//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "ChildGrafica2D.h"
#include "Clases.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "ChildWin"
#pragma resource "*.dfm"
TMDIChild3 *MDIChild3;
//---------------------------------------------------------------------------
__fastcall TMDIChild3::TMDIChild3(TComponent* Owner)
        : TMDIChild(Owner)
{
}
//---------------------------------------------------------------------------

void TMDIChild3::Grafica3D(float **F,int nr,int nc)
{
  int i,j,C=nr,R=nc,rango;
  int x0=100,y0=100;
  float Tam=1.0;
  TPoint *P;
  P=new TPoint[2*C];
  Tam=Tam*140.0/((float)nr);// Tam modify the Size of the graph
  rango=nr/2;
  for(i=-rango;i<rango;i++)
    for(j=-rango;j<rango;j++)
        Superficie[i][j]=F[i][j];
  for(j=0;j<R-1;j++)
  {
	for(i=0;i<C;i++)
	 {
		  P[i].x=x0+0.8*(-j+2.0*i)*Tam+50;
                  P[i].y=y0-F[j-rango][i-rango]*fact+Tam*j;
		  P[i+C].x=x0+0.8*(-j-1+2.0*(C-1-i))*Tam+50;
                  P[i+C].y=y0-F[j+1-rango][C-1-i-rango]*fact+Tam*(j+1);
  	 }
   for(i=0;i<(2*C)-1;i++)
    {
      Image1->Canvas->MoveTo(P[i].x,P[i].y);
      Image1->Canvas->LineTo(P[i+1].x,P[i+1].y);
    }
  }
}

void TMDIChild3::ReGrafica(void)
{
  TPoint *P;
  int rango,i,j;
  int x0=100,y0=100;
  float Tam=1.0;
  P=new TPoint[2*N];
  Tam=Tam*140.0/((float)N);// Tam modify the Size of the graph
  rango=N/2;
  //Image1->Canvas->Brush->Color = clWhite;
  Image1->Canvas->FillRect(Rect(0,0,Image1->Width,Image1->Height));
  for(j=0;j<N-1;j++)
  {
	for(i=0;i<N;i++)
	 {
		  P[i].x=x0+0.8*(-j+2.0*i)*Tam+50;
                  P[i].y=y0-Superficie[j-rango][i-rango]*fact+Tam*j;
		  P[i+N].x=x0+0.8*(-j-1+2.0*(N-1-i))*Tam+50;
                  P[i+N].y=y0-Superficie[j+1-rango][N-1-i-rango]*fact+Tam*(j+1);
  	 }
   for(i=0;i<(2*N)-1;i++)
    {
      Image1->Canvas->MoveTo(P[i].x,P[i].y);
      Image1->Canvas->LineTo(P[i+1].x,P[i+1].y);
    }
  }
}

void __fastcall TMDIChild3::Image1MouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
{
if(Button==0)
 {
  fact++;
  ReGrafica();
 }
if(Button==1)
 {
  fact--;
  ReGrafica();
 }
}
//---------------------------------------------------------------------------
void __fastcall TMDIChild3::Inicializa_Child(int factor,int rango,int Ene)
{
  Height=factor;
  Width=factor;
  Superficie=matrix(-rango,rango-1,-rango,rango-1);
  fact=3.0;
  N=Ene;
}
