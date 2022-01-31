//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include <math.h>
#include <stdlib.h>
#include <math.h>
#include "Clases.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
//---------------------------------------------------------------------------
#pragma package(smart_init)

//FUNCIONES numerical recipies

void fourn(float data[],unsigned long nn[],int ndim,int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP
/********************************************/
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
{
	int i,k;
	float p,qn,sig,un,*u;

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}
/********************************************/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
	int klo,khi,k;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0)
         ShowMessage("Error en la rutina splint, xa de entrada no vlaido");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

//FUNCIONES

void Almacena_Matriz(float **FVR,float **FReal,int rango)
{
int i,j;
for(i=-rango;i<=rango-1;i++)   //almacenar funcion pulso para graficar
  for(j=-rango;j<=rango-1;j++)
      FVR[i][j]=FReal[i][j];
}
/********************************************/
void Copia_Ms(float **FReal,float **FImag,float **Real,float **Imag,int rango)
{
int i,j;
for(i=-rango;i<=rango-1;i++)   //almacenar matrices
  for(j=-rango;j<=rango-1;j++)
      {
        FReal[i][j]=Real[i][j];
        FImag[i][j]=Imag[i][j];
      }
}
/********************************************/
void ConvR(float **f,float **g,float *gk,int nk,int n,int ir)
// convolucion de renglon ir de f con kernel 1-D
{
    int i,j,j0,j1 ;
    float s ;
	 for (i = 0 ; i < n ; i++)
    {
		s = 0 ;
      j0 = i - nk ; if (j0 < 0) j0 = 0 ;
      j1 = i + nk ; if (j1 > n-1) j1 = n-1 ;
      for (j = j0 ; j <= j1 ; j++)
         s += f[ir][j]*gk[i-j] ;
      g[ir][i] = s ;
    }
}
/****************************************************************/
void ConvC(float **f,float **g,float *gk,int nk,int n,int jc)
// convolucion de columna jc de f con kernel 1-D
{
    int i,j,j0,j1 ;
    float s ;
    for (i = 0 ; i < n ; i++)
    {
      s = 0 ;
      j0 = i - nk ; if (j0 < 0) j0 = 0 ;
      j1 = i + nk ; if (j1 > n-1) j1 = n-1 ;
      for (j = j0 ; j <= j1 ; j++)
         s += f[j][jc]*gk[i-j] ;
      g[i][jc] = s ;
	 }
}
/****************************************************************/
void Conv2D(float **g,float **f,float *gkr,float *gkc,int nk, int nr,int nc)
// convolucion 2D con kernel cuadrado
{
	 int i ;
    float **a ;
    a = matrix(0,nr-1,0,nc-1) ;
    for (i = 0 ; i < nr ; i++)
      ConvR(g,a,gkr,nk,nc,i) ;
    for (i = 0 ; i < nc ; i++)
      ConvC(a,f,gkc,nk,nr,i) ;
	 free_matrix(a,0,nr-1,0,nc-1) ;
}
/****************************************************************/
void Crea_Rectangulo(float **rectangulo,float xrec,float yrec,float num_splines,float valor)
{
int i,j,ren,col,epr,epc;
float **rectan;
ren=xrec+num_splines;
col=yrec+num_splines ;
epr=ren/2;
epc=col/2;
rectan=matrix(0,ren-1,0,col-1);
for(i=-ren/2;i<ren/2;i++)
  for(j=-col/2;j<col/2;j++)
    {
      if(i>-xrec/2 && i<xrec/2 && j>-yrec/2 && j<yrec/2)
         rectan[i+epr][j+epc]=valor;
      else
         rectan[i+epr][j+epc]=0.0;
    }
SuavizaGaussiana(rectan,rectangulo,num_splines,ren,col);
free_matrix(rectan,0,ren-1,0,col-1);
}
/********************************************/
void Datos_Archivo(FILE *APresion2D,float **Pres,int rango)
{
int i,j;
AnsiString linea;

fprintf(APresion2D,"[");
for(i=-rango;i<=rango-1;i++)
   {
     linea="";
     for(j=-rango;j<=rango-1;j++)
       {
         linea+=FloatToStr(Pres[i][j]);
         if(j!=rango-1)
           linea+=", ";
       }
     linea+=";";
     linea+='\n';
     fprintf(APresion2D,linea.c_str());
   }
fprintf(APresion2D,"]");
}
/********************************************/
float Funcion_Circulo2DS(float ra,float rb,float tex,float tey,int n_s,float ep0,float ep1,float *xa,float *ya,float *sy,float valor)  //funcion V con splines
{
float res,xs,ys,val;

xs=tex/ra;
ys=tey/rb;
val=(xs*xs)+(ys*ys);
if(val<1-ep0)
  res=valor;
else if(val>=1-ep0 && val<1+ep1 )
  splint(xa,ya,sy,n_s,val,&res);
else
  res=0.0;
return res;

}
/********************************************/
float Funcion_RectanguloS(float **rectangulo,float xrec,float yrec,int tex,int tey,float num_splines)
{
float res;
int intervalox=xrec/2,intervaloy=yrec/2;

if(fabs(tex)<intervalox+num_splines && fabs(tey)<intervaloy+num_splines)
    res=rectangulo[int(tex)+intervalox][int(tey)+intervaloy];
else
    res=0.0;
return res;
}
/********************************************/
void Funcion_Integral(complex<float> *Integral1,complex<float> *Integral2,float k,float ka,float kb,float miu,float z)
{
complex <float> K(0.0,0.0);
complex <float> H(0.0,0.0);
complex <float> Miu(0.0,0.0);
complex <float> X3(0.0,0.0);
complex <float> Ny1(0.0,0.0);
complex <float> Ny2(0.0,0.0);
complex <float> Fk(0.0,0.0);
complex <float> Dos(2.0,0.0);
complex <float> Cuatro(4.0,0.0);
K.real()=k;
Miu.real()=miu;
X3.real()=z;
H.real()=(2*k*k)-(kb*kb);
if(k>=ka)               //valores reales o imaginarios para las raices cuadras Niu
   Ny1.real()=sqrt((k*k)-(ka*ka));
else
   Ny1.imag()=sqrt((ka*ka)-(k*k));
if(k>=kb)
   Ny2.real()=sqrt((k*k)-(kb*kb));
else
   Ny2.imag()=sqrt((kb*kb)-(k*k));
Fk=(H*H)-(Cuatro*Ny1*Ny2*K*K);
*Integral1=((exp(-Ny1*X3)/(Miu*Fk))*Ny1*H);
*Integral2=-((exp(-Ny2*X3)/(Miu*Fk))*Dos*Ny1*K*K);
}
/********************************************/
void GKernel(float *gk,float sigma,int n)
// genera kernel gaussiano
{
    int i;
    float s = 0, s2 = 2.0*sigma*sigma ;
    gk[0] = 1 ;
    for (i = 1 ; i <= n ; i++)
        s += (gk[i] = exp(-(float)(i*i)/s2)) ;

    s2 = 1.0 + 2.0*s ;
    gk[0] /= s2 ;
    for (i = 1 ; i <= n ; i++)
        gk[-i] = (gk[i] /= s2) ;
}
/********************************************************/
void Genera_Splines(float *x,float *y,float *sy,int num_splines,float eps0,float eps1)
{
float yp1,ypn,uno_cero,incremen,incremen_ac=0,decremen,dif;
int i;

yp1=1e30;
ypn=1e30;
uno_cero=1;
dif=(1+eps1)-(1-eps0);  //rango para splines en X
incremen=dif/float(num_splines);     //incremento para spline en X
decremen=1/float(num_splines);    //decremento para spline en Y
for(i=1;i<=num_splines;i++,uno_cero=uno_cero-decremen)  //ciclo para inicializar tabla XY
  {
      x[i]=(1-eps0)+incremen_ac;
      incremen_ac+=incremen;
      y[i]=uno_cero;
  }
spline(x,y,num_splines,yp1,ypn,sy);   //llama a spline para crear tabla sy de resultado interpolados
}
/********************************************/
void Multiplica_Complex(float **FFR,float** FFI,float **FReal,float **FImag,float **GReal,float **GImag,int N)
{
int i,j,rango;
rango=N/2;
for(i=-rango;i<=rango-1;i++)
  for(j=-rango;j<=rango-1;j++)
    {
      FFR[i][j]=(FReal[i][j]*GReal[i][j])-(FImag[i][j]*GImag[i][j]);
      FFI[i][j]=(GReal[i][j]*FImag[i][j])+(FReal[i][j]*GImag[i][j]);
    }
}
/********************************************/
void Normas(float **F1,float **F2,int N,float *suma,float *maximo)
{
float **Norma;
int i,j,rango;
rango=N/2;
*suma=0;
*maximo=0;
Norma=matrix(-rango,rango-1,-rango,rango-1);
for(i=1;i<=N;i++)
  for(j=1;j<=N;j++)
   {
      Norma[i-rango-1][j-rango-1]=F1[i-rango-1][j-rango-1]-F2[i-rango-1][j-rango-1];
      //normas
      if(fabs(Norma[i-rango-1][j-rango-1])>*maximo)
        *maximo=fabs(Norma[i-rango-1][j-rango-1]);
      *suma+=Norma[i-rango-1][j-rango-1]*Norma[i-rango-1][j-rango-1];
   }
*suma=sqrt(*suma);
free_matrix(Norma,-rango,rango-1,-rango,rango-1);
}
/********************************************/
void Obten_Funcion_V(float **FReal,float **FImag,float delta_te,float eps0,float eps1,int num_splines,
float r,int rango,float *sy,float *x,float *y,float a,float b,float xrec,float yrec,int tipo)
{
 int i,j;
 float tei,tej,**rectangulo,valor=1.0,n_div;
 if(tipo==2)     //rectangulares
  {
   n_div=float(num_splines)*(1/float(num_splines));
   rectangulo=matrix(0,xrec+n_div,0,yrec+n_div);
  Crea_Rectangulo(rectangulo,xrec,yrec,n_div,valor);
  }
 for(i=-rango;i<=rango-1;i++)   //almacenar funcion pulso y funcion de green
   for(j=-rango;j<=rango-1;j++)
     {
       tei=i*delta_te;
       tej=j*delta_te;
       FImag[i][j]=0.0;
       if(tipo==0)     //circulares
            FReal[i][j]=Funcion_Circulo2DS(r,r,tei,tej,num_splines,eps0,eps1,x,y,sy,valor);
       if(tipo==1)     //elipticos
            FReal[i][j]=Funcion_Circulo2DS(a,b,tei,tej,num_splines,eps0,eps1,x,y,sy,valor);
       if(tipo==2)     //rectangulares
            FReal[i][j]=Funcion_RectanguloS(rectangulo,xrec,yrec,tei,tej,n_div);
     }
 if(tipo==2)     //rectangulares
   free_matrix(rectangulo,0,xrec+n_div,0,yrec+n_div);
}
/********************************************/
void Obten_TFuncion_GReen(float **GReal,float **GImag,float delta_omega,int rango,float ro,float z)
{
int i,j;
float k0,k1,k2,tek1,tek2;
complex<float> complejo(0.0,0.0);

for(i=-rango;i<=rango-1;i++)   //almacenar funcion pulso y funcion de green
  for(j=-rango;j<=rango-1;j++)
     {
        tek1=i*delta_omega;
        tek2=j*delta_omega;
        k0=sqrt((tek1*tek1)+(tek2*tek2));
        k1=k0-(delta_omega/2);
        k2=k0+(delta_omega/2);
        TDF_Green_P(&complejo,ro,z,k1,k2);
        GReal[i][j]=complejo.real();
        GImag[i][j]=complejo.imag();
     }
}
/********************************************/
void Obten_TFuncion_GReenSS(float **GReal,float **GImag,float **GReal2,float **GImag2,float delta_omega,int rango,float dens,float z,float alfa,float beta,float f,int si_rango)
{
  int i,j,longitud;
  float mu,w,k,ka,kb,tek1,tek2,h,k1,k2;
  float *Raices;
  complex<float> complejo1(0.0,0.0);
  complex<float> complejo2(0.0,0.0);
  Raices=new float[3];
  if(si_rango==1)
   w=f;
  else
   w=2*M_PI*f;
  mu=dens*beta*beta;
  ka=w/alfa;
  kb=w/beta;
  longitud=Obtener_Raices(Raices,ka,kb);
  for(i=-rango;i<=rango-1;i++)   //almacenar funcion pulso y funcion de green
    for(j=-rango;j<=rango-1;j++)
       {
        complejo1.real()=0.0;
        complejo1.imag()=0.0;
        complejo2.real()=0.0;
        complejo2.imag()=0.0;
        tek1=i*delta_omega;
        tek2=j*delta_omega;
        if((tek1*tek1)+(tek2*tek2)>0)
         {
          k=sqrt((tek1*tek1)+(tek2*tek2));
          k1=k-(delta_omega/2);
          k2=k+(delta_omega/2);
          TDF_Green_PSS(&complejo1,&complejo2,Raices,longitud,k1,k2,ka,kb,mu,z);
         }
        else  //posiscion i=0 j=0
         TDF_Green_P_S0(&complejo1,&complejo2,ka,kb,mu,z);
        GReal[i][j]=complejo1.real();
        GImag[i][j]=complejo1.imag();
        GReal2[i][j]=complejo2.real();
        GImag2[i][j]=complejo2.imag();
       }
delete Raices;
}
/********************************************/
void Obten_Presion(float **Pres,float **FFR,float **FFI,float **FReal,float **FImag,float **GReal,float **GImag,int N,float r,float a,float b,float xrec,float yrec,int transductor,int rango)
{
int i,j;
Multiplica_Complex(FFR,FFI,FReal,FImag,GReal,GImag,N);
TFourier(FFR,FFI,N,r,a,b,xrec,yrec,transductor,-1);
for(i=-rango;i<=rango-1;i++)   //presion
  for(j=-rango;j<=rango-1;j++)
      Pres[i][j]=sqrt((FFR[i][j]*FFR[i][j])+(FFI[i][j]*FFI[i][j]));
}
/********************************************/
int Obtener_Raices(float *Raices,float ka, float kb)
{
float a,b,c,d,b1,c1,chi,KTemp,R1,tol,tolerancia=0.01,funcion,derivada,longitud;
a=16*((ka*ka)-(kb*kb));
b=8*kb*kb*((3*kb*kb)-(2*ka*ka));
c=-8*pow(kb,6);
d=pow(kb,8);
R1=Obtener_Val_Ini(a,b,c,d);//valor inicial para el NR
while (tol>tolerancia)  //NR para obtener primera raiz del polinomio de grado 3
 {
   funcion=(a*R1*R1*R1)+(b*R1*R1)+(c*R1)+d;
   derivada=(3*a*R1*R1)+(2*b*R1)+c;
   KTemp=R1;
   R1=R1-(funcion/derivada);
   tol=fabs(KTemp-R1);
 }
//chicharronero
b1=b+(a*R1);
if(R1!=0)
 c1=-d/R1;
else
 c1=0.0;
if(R1<0)   //es compleja la primera raiz???
  Raices[2]=-1;
else
  Raices[2]=sqrt(R1);
chi=(b1*b1)-(4*a*c1);
if(chi<0)       //es compleja la segunda raiz???
 {
  Raices[1]=-1;
  Raices[0]=-1;
 }
else
  {
    chi=sqrt(chi);
    if(a!=0)
      R1=(-b1+chi)/(2*a);
    else
      R1=0.0;
    if(R1<0)   //es compleja la raiz???
      Raices[1]=-1;
    else
      Raices[1]=sqrt(R1);
    if(a!=0)
     R1=(-b1-chi)/(2*a);
    else
     R1=0.0;
    if(R1<0)   //es compleja la raiz???
      Raices[0]=-1;
    else
      Raices[0]=sqrt(R1);
  }
longitud=Ordena_Raices(Raices);//ordena raices
return longitud;
}
/********************************************/
float Obtener_Val_Ini(float a,float b,float c,float d)//valor inicial para el NR
{
float Val1,Val2,Val3,Res;
if(a!=0)
 {
  Val1=5*fabs(b)/fabs(a);
  Val2=sqrt(5*fabs(c)/fabs(a));
  Val3=pow(5*fabs(d)/fabs(a),1.0/3.0);
 }
else
 {
  Val1=0.0;
  Val2=0.0;
  Val3=0.0;
 }
Res=Val1;
if(Res<Val2)
 {
  Res=Val2;
  if(Res<Val3)
    Res=Val3;
 }
else if(Res<Val3)
  Res=Val3;
return (Res);
}
/********************************************/
int Ordena_Raices(float *Raices)//ordena raices
{
int i,Res=3;
float temp;
if(Raices[0]>Raices[1])
  {
    temp=Raices[1];
    Raices[1]=Raices[0];
    Raices[0]=temp;
  }
if(Raices[1]>Raices[2])
 {
   temp=Raices[2];
   Raices[2]=Raices[1];
   Raices[1]=temp;
 }
if(Raices[0]>Raices[1])
 {
    temp=Raices[1];
    Raices[1]=Raices[0];
    Raices[0]=temp;
 }
for(i=0;i<3;i++)
 {
   if(Raices[i]<0)
     Res--;
 }
return Res;  //regresa el numero de raices positivas, ordenadas de der, a izq en el vector raices
}
/********************************************/
void Parametros(int N,float r,float a,float b,float xrec, float yrec,float *A,float *delta_te,float *delta_omega,int *rango,int tipo)
{
switch(tipo)
 {
   case 0:
        *A=4*r;
        break;
   case 1:
        *A=4*sqrt(a*b);
        break;
   case 2:
        *A=4*sqrt(xrec*yrec);
        break;
 }
*delta_te=*A/float(N);
*delta_omega=2*M_PI/(*A);
*rango=N/2;
}
/********************************************/
int Raices_Dentro(float *Raices_in,float *Raices,int longitud,float k1,float k2)
{
 int i,num=0;

 for(i=2;i>=3-longitud;i--)
  {
    if(Raices[i]<=k2 && Raices[i]>=k1)
     {
      Raices_in[2-i]=Raices[i];
      num++;
     }
  }
return num;
}
/****************************************************************/
void SuavizaGaussiana(float **dat,float **res,float sigma,int nr,int nc)
//*  Suaviza con kernel gaussiano con parametro sigma
{
    float *gk ;
    int nk = 3.0*sigma ;
    if(nk)
    {
      gk = vector(-nk,nk) ;
		GKernel(gk,sigma,nk) ;
      Conv2D(dat,res,gk,gk,nk,nr,nc) ;
		free_vector(gk,-nk,nk) ;
	 }
    else
    {
      for (int i = 0 ; i < nr ; i++)
        for (int j = 0 ; j < nc ; j++)
           res[i][j] = dat[i][j];
	 }
}
/********************************************/
void TDF_Green_P(complex<float> *complejo, float ro,float d,float k1,float k2)
{
complex <float> I(0.0,1.0);
complex <float> D(0.0,0.0);
complex <float> K(0.0,0.0);
complex <float> Ge(0.0,0.0);
complex <float> Ro(0.0,0.0);
complex <float> W(0.0,0.0);
complex <float> Dos(2.0,0.0);
float f=3.5,c=1.4883,w,k,k12,k22,ka2,dif2;
w=2*M_PI*f;
k=w/c;
D.real()=d;
K.real()=k;
W.real()=w;
Ro.real()=ro;
k12=k1*k1;
k22=k2*k2;
ka2=k*k;
dif2=k22-k12;
   //3 casos para d=0
if(dif2!=0)
 {
    if(d==0)
     {
       if(k<k1)
         {
          Ge.real()=2*((sqrt(k22-ka2)-sqrt(k12-ka2))/dif2);
          Ge.imag()=0.0;
         }
       if(k1<=k && k<=k2)
         {
          Ge.real()=2*(sqrt(k22-ka2)/dif2);
          Ge.imag()=2*(sqrt(ka2-k12)/dif2);
         }
       if(k>k2)
         {
          Ge.imag()=2*((sqrt(ka2-k12)-sqrt(ka2-k22))/dif2);
          Ge.real()=0.0;
         }
     }
   //3 casos para d>0
    if(d>0)
     {
       if(k<k1)
         {
          Ge.real()=2*((exp(-sqrt(k12-ka2)*d)-exp(-sqrt(k22-ka2)*d))/(dif2*d));
          Ge.imag()=0.0;
         }
       if(k1<=k && k<=k2)
         {
          Ge.real()=2*((cos(sqrt(ka2-k12)*d)-exp(-sqrt(k22-ka2)*d))/(dif2*d));
          Ge.imag()=2*(sin(sqrt(ka2-k12)*d)/(dif2*d));
         }
       if(k>k2)
         {
          Ge.real()=2*((cos(sqrt(ka2-k12)*d)-cos(sqrt(ka2-k22)*d))/(dif2*d));
          Ge.imag()=2*((sin(sqrt(ka2-k12)*d)-sin(sqrt(ka2-k22)*d))/(dif2*d));
         }
     }
 }
else
  Ge=Dos*I*exp(I*K*D)/K;
*complejo=I*W*Ro*Ge;
}
/********************************************/
void TDF_Green_P_S0(complex<float> *complejo1,complex<float> *complejo2,float ka,float kb,float mu,float z)
{
complex <float> Ny1(0.0,0.0);
complex <float> X3(0.0,0.0);
complex <float> Miu(0.0,0.0);
complex <float> H(0.0,0.0);
complex <float> KBeta4(0.0,0.0);
X3.real()=z;
Miu.real()=mu;
H.real()=-(kb*kb);
Ny1.real()=0.0;
Ny1.imag()=ka;
KBeta4.real()=pow(kb,4);
*complejo1=((exp(-Ny1*X3)/(Miu*KBeta4))*Ny1*H);
*complejo2=0.0;
}
/********************************************/
void TDF_Green_PSS(complex<float> *complejo1,complex<float> *complejo2,float *Raices,int longitud,float k1,float k2,float ka,float kb,float miu,float z)
{
complex <float> Paso(0.0,0.0);
complex <float> Integral1(0.0,0.0);
complex <float> Integral2(0.0,0.0);
complex <float> UnoCinco(1.5,0.0);
complex <float> PuntoCinco(0.5,0.0);
float paso,*Raices_in,h;
int i,num_raices,N=5;
//integrales
Raices_in=new float[3];
num_raices=Raices_Dentro(Raices_in,Raices,longitud,k1,k2);
if(num_raices>=1)
 {
   paso=(Raices_in[0]-k1)/N;
   Paso.real()=paso;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[0],ka,kb,miu,z);
   *complejo1=UnoCinco*Paso*Integral1;
   *complejo2=UnoCinco*Paso*Integral2;
   Funcion_Integral(&Integral1,&Integral2,k1,ka,kb,miu,z);
   *complejo1+=PuntoCinco*Paso*Integral1;
   *complejo2+=PuntoCinco*Paso*Integral2;
   for(i=1;i<N-1;i++)
     {
       Funcion_Integral(&Integral1,&Integral2,k1+(i*paso),ka,kb,miu,z);
       *complejo1+=paso*Integral1;
       *complejo2+=paso*Integral2;
     }
 }
if(num_raices>=2)
 {
   paso=(Raices_in[1]-Raices_in[0])/N;
   Paso.real()=paso;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[0],ka,kb,miu,z);
   *complejo1+=UnoCinco*Paso*Integral1;
   *complejo2+=UnoCinco*Paso*Integral2;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[1],ka,kb,miu,z);
   *complejo1+=UnoCinco*Paso*Integral1;
   *complejo2+=UnoCinco*Paso*Integral2;
   for(i=1;i<N-1;i++)
     {
       Funcion_Integral(&Integral1,&Integral2,Raices_in[0]+(i*paso),ka,kb,miu,z);
       *complejo1+=paso*Integral1;
       *complejo2+=paso*Integral2;
     }
 }
if(num_raices>=3)
 {
   paso=(Raices_in[2]-Raices_in[1])/N;
   Paso.real()=paso;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[1],ka,kb,miu,z);
   *complejo1+=UnoCinco*Paso*Integral1;
   *complejo2+=UnoCinco*Paso*Integral2;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[2],ka,kb,miu,z);
   *complejo1+=UnoCinco*Paso*Integral1;
   *complejo2+=UnoCinco*Paso*Integral2;
   for(i=1;i<N-1;i++)
     {
       Funcion_Integral(&Integral1,&Integral2,Raices_in[1]+(i*paso),ka,kb,miu,z);
       *complejo1+=paso*Integral1;
       *complejo2+=paso*Integral2;
     }
 }
if(num_raices>0)
 {
   paso=(k2-Raices_in[num_raices-1])/N;
   Paso.real()=paso;
   Funcion_Integral(&Integral1,&Integral2,k2,ka,kb,miu,z);
   *complejo1+=PuntoCinco*Paso*Integral1;
   *complejo2+=PuntoCinco*Paso*Integral2;
   Funcion_Integral(&Integral1,&Integral2,Raices_in[num_raices-1],ka,kb,miu,z);
   *complejo1+=UnoCinco*Paso*Integral1;
   *complejo2+=UnoCinco*Paso*Integral2;
   for(i=1;i<N-1;i++)
     {
       Funcion_Integral(&Integral1,&Integral2,Raices_in[num_raices-1]+(i*paso),ka,kb,miu,z);
       *complejo1+=paso*Integral1;
       *complejo2+=paso*Integral2;
     }
 }
if(num_raices==0)
 {
   paso=(k2-k1)/N;
   Paso.real()=paso;
   Funcion_Integral(&Integral1,&Integral2,k2,ka,kb,miu,z);
   *complejo1=PuntoCinco*Paso*Integral1;
   *complejo2=PuntoCinco*Paso*Integral2;
   Funcion_Integral(&Integral1,&Integral2,k1,ka,kb,miu,z);
   *complejo1+=PuntoCinco*Paso*Integral1;
   *complejo2+=PuntoCinco*Paso*Integral2;
   for(i=1;i<N-1;i++)
     {
       Funcion_Integral(&Integral1,&Integral2,k1+(i*paso),ka,kb,miu,z);
       *complejo1+=paso*Integral1;
       *complejo2+=paso*Integral2;
     }
 }
delete Raices_in;
}
/********************************************/
void TFourier(float **Real,float **Imag,int N,float r,float a,float b,float x, float y,int transductor,int tipo)
{
unsigned long *nlong;
float *hn;
float A,delta_te,delta_omega;
int i,j,rango;

hn=vector(1,(2*N*N));
nlong=new unsigned long[3];  //dimensiones 3
nlong[1]=N;
nlong[2]=N;
Parametros(N,r,a,b,x,y,&A,&delta_te,&delta_omega,&rango,transductor);
if(tipo==1)     //directa
 {
    for(i=1;i<=N;i++)     //arregla datos para el vector de numerical
      for(j=1;j<=N;j++)
        {
          hn[((i-1)*2*N)+((2*j)-1)]=pow(-1,(i+j))*Real[i-rango-1][j-rango-1]*delta_te*delta_te;
          hn[((i-1)*2*N)+(2*j)]=pow(-1,(i+j))*Imag[i-rango-1][j-rango-1]*delta_te*delta_te;
        }
    fourn(hn,nlong,2,-1);   //FOURIER F Directa
    for(i=1;i<=N;i++)   // Obtiene Fkl
      for(j=1;j<=N;j++)
        {
          Real[i-rango-1][j-rango-1]=pow(-1,(i+j))*hn[((i-1)*2*N)+((2*j)-1)];
          Imag[i-rango-1][j-rango-1]=pow(-1,(i+j))*hn[((i-1)*2*N)+(2*j)];
        }
 }
if(tipo==-1)  //inversa
 {

    for(i=1;i<=N;i++)    //Para obtener inversa, arreglar Fkl
      for(j=1;j<=N;j++)
        {
           hn[((i-1)*2*N)+((2*j)-1)]=pow(-1,(i+j))*Real[i-rango-1][j-rango-1]/(A*A);
           hn[((i-1)*2*N)+(2*j)]=pow(-1,(i+j))*Imag[i-rango-1][j-rango-1]/(A*A);
        }
    fourn(hn,nlong,2,1);  //obtiene FOURIER F Inversa de hn en espacio de fourier
    for(i=1;i<=N;i++)
      for(j=1;j<=N;j++)     //Obtiene fmn
        {
          Real[i-rango-1][j-rango-1]=pow(-1,(i+j))*hn[((i-1)*2*N)+((2*j)-1)];
          Imag[i-rango-1][j-rango-1]=pow(-1,(i+j))*hn[((i-1)*2*N)+(2*j)];
        }
 }
delete nlong;
free_vector(hn,1,(2*N*N));
}

/********************************************/
//funciones de memoria

float *vector(int i0,int i1)
/* allocates memory for vector with indices from i0 to i1 */
{
	 float *v ;
	if ((v = (float *)calloc(i1-i0+1,sizeof (float))) == NULL)
        ShowMessage("error de memoria") ;
	 return v-i0 ;
}
/***********************************************/
void free_vector(float *v,int i0, int i1)
/* frees memory for vector v */
{
	free(v + i0) ;
}
void free_vector(float *v,int i0)
/* frees memory for vector v */
{
	free(v + i0) ;
}
float **matrix(int nrl,int nrh,int ncl,int nch)
{
    int i;
    float **m;
    m= new float *[nrh-nrl+1];
    if (!m) ShowMessage("fallo la asignaci¢n de memoria en el paso 1 par la matriz()");
    m -= nrl;
    for (i=nrl;i<=nrh;i++)
    {
      m[i]=new float [nch-ncl+1];
              if (!m[i]) ShowMessage("fallo la asignaci¢n de memoria en el paso 2 par la matriz()");
                m[i] -= ncl;
    }
    return m;
}
/********************************************/
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch)
{
    int i;
    for(i=nrh;i>=nrl;i--)
        delete[] (m[i]+ncl);
    delete[] (m+nrl);
}

