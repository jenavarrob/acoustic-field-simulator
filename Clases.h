//---------------------------------------------------------------------------
#ifndef ClasesH
#define ClasesH
#include <complex.h>
//---------------------------------------------------------------------------

//funciones numerical recipies
void fourn(float data[],unsigned long nn[],int ndim,int isign);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
//funciones
void Almacena_Matriz(float **FVR,float **FReal,int rango);
void Copia_Ms(float **FReal,float **FImag,float **Real,float **Imag,int rango);
void ConvR(float **f,float **g,float *gk,int nk,int n,int ir);
void ConvC(float **f,float **g,float *gk,int nk,int n,int jc);
void Conv2D(float **g,float **f,float *gkr,float *gkc,int nk, int nr,int nc);
void Crea_Rectangulo(float **rectangulo,int xrec,int yrec,float num_splines,float valor);
void Datos_Archivo(FILE *APresion2D,float **Pres,int rango);
float Funcion_Circulo2DS(float ra,float rb,float tex,float tey,int n_s,float ep0,float ep1,float *xa,float *ya,float *sy,float valor); //funcion V con splines
float Funcion_RectanguloS(float **rectangulo,float xrec,float yrec,int tex,int tey,float num_splines);
void Funcion_Integral(complex<float> *Integral1,complex<float> *Integral2,float k,float ka,float kb,float miu,float z);
void GKernel(float *gk,float sigma,int n);
void Genera_Splines(float *x,float *y,float *sy,int num_splines,float eps0,float eps1);
void Multiplica_Complex(float **FFR,float** FFI,float **FReal,float **FImag,float **GReal,float **GImag,int N);
void Normas(float **F1,float **F2,int N,float *suma,float *maximo); //obtiene norma de dos matrices
void Obten_Funcion_V(float **FReal,float **FImag,float delta_te,float eps0,float eps1,int num_splines,
float r,int rango,float *sy,float *x,float *y,float a,float b,float xrec,float yrec,int tipo);
void Obten_TFuncion_GReen(float **Greal,float **GImag,float delta_omega,int rango,float ro,float z);
void Obten_TFuncion_GReenSS(float **GReal,float **GImag,float **GReal2,float **GImag2,float delta_omega,int rango,float dens,float z,float alfa,float beta,float f,int si_rango);
void Obten_Presion(float **Pres,float **FFR,float **FFI,float **FReal,float **FImag,float **GReal,float **GImag,int N,float r,float a,float b,float xrec,float yrec,int transductor,int rango);
int Obtener_Raices(float *Raices,float ka, float kb);
float Obtener_Val_Ini(float a,float b,float c,float d);//valor inicial para el NR
int Ordena_Raices(float *Raices);//ordena raices
void Parametros(int N,float r,float a, float b,float xrec, float yrec,float *A,float *delta_te,float *delta_omega,int *rango,int tipo); //obtiene parametros de fourier
int Raices_Dentro(float *Raices_in,float *Raices,int longitud,float k1,float k2);
void SuavizaGaussiana(float **dat,float **res,float sigma,int nr,int nc);
void TDF_Green_P(complex<float> *complejo,float ro,float z,float tei,float tej);//TDF de funcion de green sin singularidad conversion a polar
void TDF_Green_P_S0(complex<float> *complejo1,complex<float> *complejo2,float ka,float kb,float mu,float z);
void TDF_Green_PSS(complex<float> *complejo1,complex<float> *complejo2,float *Raices,int longitud,float k1,float k2,float ka,float kb,float miu,float z);
void TFourier(float **Real,float **Imag,int N,float r,float a,float b,float x, float y,int transductor,int tipo);  //TDF con modificacion deseada para acusticos
//funciones memoria
float *vector(int i0,int i1);
void free_vector(float *v,int i0,int i1);
void free_vector(float *v,int i0);
float **matrix(int nrl,int nrh,int ncl,int nch);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);

#endif
