#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include "lodepng.h"
#include <fftw3.h>
#include <png.h>

/*
 * Compile line: 
 * gcc -O3 -pipe dtensor_c.c uc_alone.c -o uc.elf -lm -lfftw3f -lfftw3 -lfftw3f_threads -lfftw3_threads -lpng
 * */

unsigned int NUMBER_OF_THREADS=8;


double fnxx (double x,double y,double z,double dx,double dy,double dz);
double fnxy (double x,double y,double z,double dx,double dy,double dz);
double Dz   (double a,double b,double c);


void run (int steps);
void analytic(int steps);
void hc_problem_1();
void init();
void finish();
void heff();
void demag_tensor();
void LLG(float *k);
void onestep();
void pngsave();
void savepng32 (int x,int y,void *ptr,char *filename);

const double pi = 3.14159265358979323846264338327950288419716939937511;
const double mu0= 1.25663706143591729538505735331180115367886775975004e-6 ;
const double gamma0 = 2.211e5 ;

fftwf_plan dmgf,dmgb;
fftwf_complex *cm,*ch,*nn;
float *rm,*tm,*rh,*bm,*k1,*k2,*k3,*k4;

int NX,NY,KX,KY,TNX;
unsigned int totalsteps;
double kvolume;
float efactor;
float delta,h,gammap,lambda,alfa;
float Ms,A,Ku,Dind,tau,ttime;
float fMs;
float ftau,htau,ttau,stau;
float Edmg,Eex,Edmi,Epma;
float Fpma,Fex,Fdmi;
float Kex;
float ukx,uky,ukz;
float ukxx,ukyy,ukzz,ukxy,ukyz,ukzx;
int tiempos[100]; // Para profiling interno con clock();
unsigned int verbose,debugging,nodemag,nohex;

int main()
{
  verbose = 0;
  debugging = 0;
  nodemag = 0;
  nohex = 0 ;
  hc_problem_1();
  pngsave();
  run(240000);
  finish();
  return 0;
}


void onestep(){
  static double ptime1=0.;
  static double ptime2=0.;
  static int k=0;
  ptime1+=tau;
  ptime2+=tau;
  //pngsave();
  if (5e-12 < ptime1) {fprintf(stderr,"P");ptime1 -= 5e-12; pngsave();}
  if (5e-13 < ptime2) {k++; if (k<=9) fprintf(stderr,"."); else k=0 ;ptime2 -= 5e-13;}
  return;
}

void pngsave(){
  static unsigned int counter=0;
  char name[80];
  unsigned error;
  sprintf (name,"m%06d.png",counter);
  if (verbose) printf ("Guardando %s \n",name);
  unsigned char *bitmap;
  bitmap = malloc (4*NX*NY);
  float vx,vy,vz;
  float ur,ug,ub;
  float h, s, l; //this function works with floats between 0 and 1
  float temp1, temp2, tempr, tempg, tempb;
  float fi;
  
  for (int y=0;y<NY;y++)  for (int x=0;x<NX;x++)
  {
    vx=rm[3*(y*NX+x)+0]/fMs;
    vy=rm[3*(y*NX+x)+1]/fMs;
    vz=rm[3*(y*NX+x)+2]/fMs;
    //if (debugging) vx=vy=sqrtf(0.5f);
    //if (debugging) vz=0.f;
    if (vz>1.f)  vz=1.f;
    if (vz<-1.f) vz=-1.f;
    fi = atan2f (-vy,vx);
    if (fi<0.f) fi+=2.f*pi;
    fi = fi/(2.f*pi);
  
    h = 1.0f-fi;
    s = 1.0f;
    l = 0.5f * (vz + 1.0f);

    //set the temporary values
    if(l < 0.5f) temp2 = l * (1.f + s);
    else temp2 = (l + s) - (l * s);
    temp1 = 2.f * l - temp2;
    tempr=h + 1.0f / 3.0f;
    if(tempr > 1.0f) tempr--;
    tempg=h;
    tempb=h-1.0f / 3.0f;
    if(tempb < 0.0f) tempb++;
  

    //red
    if(tempr < 1.0f / 6.0f) ur = temp1 + (temp2 - temp1) * 6.0f * tempr;
    else if(tempr < 0.5f) ur = temp2;
    else if(tempr < 2.0f / 3.0f) ur = temp1 + (temp2 - temp1) * ((2.0f / 3.0f) - tempr) * 6.0f;
    else ur = temp1;
  
    //green
    if(tempg < 1.0 / 6.0) ug = temp1 + (temp2 - temp1) * 6.0 * tempg;
    else if(tempg < 0.5) ug=temp2;
    else if(tempg < 2.0 / 3.0) ug = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempg) * 6.0;
    else ug = temp1;
  
    //blue
    if(tempb < 1.0 / 6.0) ub = temp1 + (temp2 - temp1) * 6.0 * tempb;
    else if(tempb < 0.5) ub = temp2;
    else if(tempb < 2.0 / 3.0) ub = temp1 + (temp2 - temp1) * ((2.0 / 3.0) - tempb) * 6.0;
    else ub = temp1; 
    // 0r 1g 2b 3a
    bitmap[4*(y*NX+x)+0]=(unsigned char) (255.f*ur);
    bitmap[4*(y*NX+x)+1]=(unsigned char) (255.f*ug);
    bitmap[4*(y*NX+x)+2]=(unsigned char) (255.f*ub);
    bitmap[4*(y*NX+x)+3]=0xff;
  }
  savepng32 (NX,NY,bitmap,name);
  //error = lodepng_encode32_file(name,bitmap,NX,NY);
  //if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  
  free(bitmap);
  counter++;
  return;
}

void analytic(int steps){
  static float vx,vy,vz,wx,wy,wz,omega,R,mx,my,mz,omegai,norma;
  double omegasum;
  for (;steps!=0;steps--)
  {
    omegasum=0;
    heff();
    for (int x=0;x<NX*NY;x++)
    {
      mx = rm[3*x+0];                 my = rm[3*x+1];                 mz = rm[3*x+2];
      vx = gammap * rh[3*x+0];        vy = gammap * rh[3*x+1];        vz = gammap * rh[3*x+2];
      wx = vx + alfa*(vz*my-vy*mz);   wy = vy + alfa*(vx*mz-vz*mx);   wz = vz + alfa*(vy*mx-vx*my);
      vx = wz*my-wy*mz;               vy = wx*mz-wz*mx;               vz = wy*mx-wx*my;
      
      omega = sqrtf(wx*wx+wy*wy+wz*wz);               R = mx*wx+my*wy+mz*wz;
      omegasum += omega;
      if (omega>1e-12) if (R>fMs)
      {
        omegai = 1.f/omegai;
        
        rm[3*x+0] = omegai*(-R*wx+(R*wx-mx*omega*omega)*cosf(omega*ftau)-vx*omega*sinf(omega*ftau));
        rm[3*x+1] = omegai*(-R*wy+(R*wy-my*omega*omega)*cosf(omega*ftau)-vy*omega*sinf(omega*ftau));
        rm[3*x+2] = omegai*(-R*wz+(R*wz-mz*omega*omega)*cosf(omega*ftau)-vz*omega*sinf(omega*ftau));
        norma = fMs/sqrtf((rm[3*x]*rm[3*x])+(rm[3*x+1]*rm[3*x+1])+(rm[3*x+2]*rm[3*x+2]));
        rm[3*x + 0] *= norma;     rm[3*x + 1] *= norma;       rm[3*x + 2] *= norma;  
      }
    }
    printf ("Pass %d , omega = %lf (%lg)\n",totalsteps,omegasum/(NX*NY),omegasum/(NX*NY));
    totalsteps++;                                                                  
    onestep();                                                                     
    ttime += tau;
  }
  return;
}

void run(int steps){

for (;steps!=0;steps--)
{
  memcpy(bm,rm,NX*NY*3*sizeof(float));
  
  heff();  LLG(k1);
  for (int x=0;x<3*NX*NY;x++) rm[x]+=htau*k1[x];
  heff();  LLG(k2);
  memcpy(rm,bm,NX*NY*3*sizeof(float));
  for (int x=0;x<3*NX*NY;x++) rm[x]+=htau*k2[x];
  heff();  LLG(k3);
  memcpy(rm,bm,NX*NY*3*sizeof(float));
  for (int x=0;x<3*NX*NY;x++) rm[x]+=htau*k3[x];
  heff();  LLG(k4);
  memcpy(rm,bm,NX*NY*3*sizeof(float));
  
  for (int x=0;x<3*NX*NY;x++) rm[x]+= stau*(k1[x]+k4[x]) + ttau*(k2[x]+k3[x]);
  
  
  float norma;
  for (int x=0;x<NX*NY;x++) 
  {
    norma = fMs/sqrtf((rm[3*x]*rm[3*x])+(rm[3*x+1]*rm[3*x+1])+(rm[3*x+2]*rm[3*x+2]));
    rm[3*x + 0] *= norma;     rm[3*x + 1] *= norma;       rm[3*x + 2] *= norma;  
  }
                                                                             
                                                                                 
  totalsteps++;                                                                  
  onestep();                                                                     
  ttime += tau;                                                                  
}                                                                                
return;                                                                          
}                                                                                
void LLG(float *k){                                                              
  float mhx,mhy,mhz;                                                             
  float mmhx,mmhy,mmhz;                                                          
  for (int x=0;x<NY*NX;x++)                                                         
  {                                                                              
    mhx  = (rm[3*x+1]*rh[3*x+2]) - (rm[3*x+2]*rh[3*x+1]);                            
    mhy  = (rm[3*x+2]*rh[3*x+0]) - (rm[3*x+0]*rh[3*x+2]);                            
    mhz  = (rm[3*x+0]*rh[3*x+1]) - (rm[3*x+1]*rh[3*x+0]);                            
                                                                                
    mmhx = (rm[3*x+1]*mhz) - (rm[3*x+2]*mhy);                                      
    mmhy = (rm[3*x+2]*mhx) - (rm[3*x+0]*mhz);                                      
    mmhz = (rm[3*x+0]*mhy) - (rm[3*x+1]*mhx);                                      
                                                                                 
    k[3*x+0] = gammap*mhx + lambda*mmhx;                                         
    k[3*x+1] = gammap*mhy + lambda*mmhy;                                         
    k[3*x+2] = gammap*mhz + lambda*mmhz;                                         
    }                                                                            
  return;                                                                        
}                                                                                
void hc_problem_1 (){                                                            
  NX = 256;                                                                      
  NY = 64;                                                                      
                                                                                 
  delta = 1e-9;                                                                  
  h     = 1.;                                                                    
  alfa = 0.1;                                                                      
  Ms = 800e3;                                                                    
  A = 15e-12;                                                                    
  Ku = 700e3;                                                                    
  Dind = 1.8e-3;                                                                 
  ukx=uky=0.; ukz=1.;                                                            
  tau  = 1e-14;                                                                  
  init();                                                                        
  float tmp,tmp2;                                                                     
  tmp = fMs*sqrtf(0.333f);    
  for (int x=0;x<3*NX*NY;x++) rm[x]=tmp;                                                   
  /*                                                                             
  for (int y=0;y<NY;y++)  for(int x=0;x<NX;x++)                                  
  {
    rm[3*(y*NX+x) +0] = 0;
    rm[3*(y*NX+x) +1] = -fMs;
    rm[3*(y*NX+x) +2] = 0;
    if (x<=y) if (x<=(NY-y))  
    { 
      rm[3*(y*NX+x) +0] = 0;
      rm[3*(y*NX+x) +1] = fMs;
      rm[3*(y*NX+x) +2] = 0; 
    }
    if (x<=y) if (y>=(NX-x))  
    { 
      rm[3*(y*NX+x) +0] = -fMs;
      rm[3*(y*NX+x) +1] = 0;
      rm[3*(y*NX+x) +2] = 0; 
    }
    if (y<=x) if (x<=(NY-y))  
    { 
      rm[3*(y*NX+x) +0] = fMs;
      rm[3*(y*NX+x) +1] = 0;
      rm[3*(y*NX+x) +2] = 0; 
    }
  }
  */
  
  /*
  for (int y=0;y<NY;y++)  for(int x=0;x<NX;x++) 
  {
    //tmp = atan2(x-(NX/2),y-(NY/2));
    tmp2 = sqrtf(powf(x-(NX/2),2.f)+powf(y-(NY/2),2.f)+NX/4);
    rm[3*(y*NX+x) +0] =  fMs*(x-(NX/2))/tmp2;
    rm[3*(y*NX+x) +1] = -fMs*(y-(NY/2))/tmp2;
    rm[3*(y*NX+x) +2] =  NX*fMs/tmp2/4;
  }
  */
  return ;
}
void init(){
  fftwf_init_threads();
  fftwf_plan_with_nthreads(NUMBER_OF_THREADS);
  
  KX = 2*NX; 
  KY = 2*NY;
  TNX = 3*NX;
  totalsteps =0;
  
  efactor = (float) NX * (float) NY * mu0 * (-0.5f) * delta * delta * delta * h; 
  
  kvolume = 1. / ((double) (KX*KY));
  
  gammap = -gamma0/(1.f +(alfa*alfa));
  lambda = gammap *alfa / Ms;
    
  ftau = (float) tau;
  htau = ftau / 2.f;
  ttau = ftau / 3.f;
  stau = ttau / 2.f; 
  
  fMs  = (float) Ms;   
  Fpma = 2*Ku/(mu0*Ms*Ms);
  Fex  = 2*A/(mu0*Ms*Ms);
  Fdmi = 2*Dind/(mu0*Ms*Ms);
  Kex  = Fex * powf(delta,-2.f);
  
  ukxx=Fpma*ukx*ukx;  ukyy=Fpma*uky*uky;  ukzz=Fpma*ukz*ukz;
  ukxy=Fpma*ukx*uky;  ukyz=Fpma*uky*ukz;  ukzx=Fpma*ukz*ukx;
  
  // Reserva de memoria
  nn = (fftwf_complex*) fftwf_malloc (6*KX*KY*sizeof(fftwf_complex));
  demag_tensor(); // Esto requiere mucha memoria mientras se ejecuta, así que se hace ahora.
  
  cm = (fftwf_complex*) fftwf_malloc (3*KX*KY*sizeof(fftwf_complex));   // Array for Hdmg
  ch = (fftwf_complex*) fftwf_malloc (3*KX*KY*sizeof(fftwf_complex));   // Array for Hdmg
  rm = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Magnetization
  rh = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Effective-field
  bm = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Backup Magnetization
  k1 = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Runge-Kutta k1
  k2 = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Runge-Kutta k2
  k3 = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Runge-Kutta k3
  k4 = (float*)         fftwf_malloc (3*NX*NY*sizeof(float));           // Runge-Kutta k4
  
  int n[2]={KY,KX};
  
  dmgf = fftwf_plan_many_dft(2,n,3,ch,n,3,1,cm,n,3,1,FFTW_FORWARD,FFTW_MEASURE|FFTW_DESTROY_INPUT);
  dmgb = fftwf_plan_many_dft(2,n,3,ch,n,3,1,cm,n,3,1,FFTW_BACKWARD,FFTW_MEASURE|FFTW_DESTROY_INPUT);
  return; 
}
void finish (){
  fftwf_free(nn);
  fftwf_free(cm);
  fftwf_free(ch);
  fftwf_free(tm);
  fftwf_free(rh);
  
  fftwf_free(bm);
  fftwf_free(k1);
  fftwf_free(k2);
  fftwf_free(k3);
  fftwf_free(k4);
  
  fftwf_destroy_plan(dmgf);
  fftwf_destroy_plan(dmgb);
  fftwf_cleanup_threads();
  //fftw_cleanup();
}

void heff(){  
  float *i,*o,*t;
  // Step 1 : Copy Real-m to Complex-H
  memset(ch,0,sizeof(fftwf_complex)*KX*KY*3);
  memset(cm,0,sizeof(fftwf_complex)*KX*KY*3);
  for (int y=0;y<NY;y++)
  {
    i = &rm[y*NX*3];
    o = (float*) ch ; o += KX*6*y;
    for (int x=0;x<3*NX;x++)
    o[2*x] = i[x];
  }
  
  // Step 2 : Perform a Forward FFT from H to M
  fftwf_execute(dmgf);
  // Step 3 : Apply Dmg tensor to M and store in H
    i = (float*) cm;
    o = (float*) ch;
    t = (float*) nn;
    for (int x=0;x<KX*KY;x++)
    {
      o[6*x + 0] = (i[6*x+0]*t[12*x+ 0]) - (i[6*x+1]*t[12*x+ 1]) + (i[6*x+2]*t[12*x+ 6]) - (i[6*x+3]*t[12*x+ 7]) + (i[6*x+4]*t[12*x+10]) - (i[6*x+5]*t[12*x+11]);
      o[6*x + 1] = (i[6*x+0]*t[12*x+ 1]) + (i[6*x+1]*t[12*x+ 0]) + (i[6*x+2]*t[12*x+ 7]) + (i[6*x+3]*t[12*x+ 6]) + (i[6*x+4]*t[12*x+11]) + (i[6*x+5]*t[12*x+10]);
      o[6*x + 2] = (i[6*x+0]*t[12*x+ 6]) - (i[6*x+1]*t[12*x+ 7]) + (i[6*x+2]*t[12*x+ 2]) - (i[6*x+3]*t[12*x+ 3]) + (i[6*x+4]*t[12*x+ 8]) - (i[6*x+5]*t[12*x+ 9]);
      o[6*x + 3] = (i[6*x+0]*t[12*x+ 7]) + (i[6*x+1]*t[12*x+ 6]) + (i[6*x+2]*t[12*x+ 3]) + (i[6*x+3]*t[12*x+ 2]) + (i[6*x+4]*t[12*x+ 9]) + (i[6*x+5]*t[12*x+ 8]);
      o[6*x + 4] = (i[6*x+0]*t[12*x+10]) - (i[6*x+1]*t[12*x+11]) + (i[6*x+2]*t[12*x+ 8]) - (i[6*x+3]*t[12*x+ 9]) + (i[6*x+4]*t[12*x+ 4]) - (i[6*x+5]*t[12*x+ 5]);
      o[6*x + 5] = (i[6*x+0]*t[12*x+11]) + (i[6*x+1]*t[12*x+10]) + (i[6*x+2]*t[12*x+ 9]) + (i[6*x+3]*t[12*x+ 8]) + (i[6*x+4]*t[12*x+ 5]) + (i[6*x+5]*t[12*x+ 4]);
    }
      
  fftwf_execute (dmgb);
  /* Michele Voto :  Normalizar el tensor y no el campo por la linealidad de la FFT */
  for (int y=0;y<NY;y++)
  {
    o = &rh[y*NX*3]; i = (float*) &cm[y*KX*3][0];
    for (int x=0;x<3*NX;x++) o[x] = i[2*x];
  }
  return;
}

void demag_tensor (){
  fftw_init_threads();
  fftw_plan_with_nthreads(NUMBER_OF_THREADS);
  
  fftw_complex *tmp,*row;
  fftw_plan tensor;
  
  int n[2]={KY,KX};
  tmp = fftw_malloc(6*KX*KY*sizeof(fftw_complex));
  tensor = fftw_plan_many_dft(2,n,6,tmp,n,6,1,tmp,n,6,1,FFTW_FORWARD,FFTW_ESTIMATE);
  
  for (int x=0;x<6*KX*KY;x++) tmp[x][0]=0.;
  for (int x=0;x<6*KX*KY;x++) tmp[x][1]=0.;
  fprintf (stderr,"Calculating Demag Tensor .");
  for (int y=0;y<KY;y++) 
  {
    row = &tmp[y*KX*6];
    for (int x=0;x<KX;x++)
    {
      double cx,cy;
      if (y>=NY) cy = (double) (y-KY) ; else cy = (double) y;
      if (x>=NX) cx = (double) (x-KX) ; else cx = (double) x;
          
      row[6*x + 0][0] = fnxx(cx,cy,0.,1.,1.,h ); // xx
      row[6*x + 1][0] = fnxx(cy,cx,0.,1.,1.,h ); // yy
      row[6*x + 2][0] = fnxx(0.,cy,cx,h ,1.,1.); // zz
      row[6*x + 3][0] = fnxy(cx,cy,0.,1.,1.,h ); // xy
      row[6*x + 4][0] = fnxy(cy,0.,cx,1.,h ,1.); // yz 
      row[6*x + 5][0] = fnxy(cx,0.,cy,1.,h ,1.); // zx 
      
    }
  }
  
  tmp[0][0] = Dz(h ,1.,1.) -2*Kex; 
  tmp[1][0] = Dz(1.,h ,1.) -2*Kex;
  tmp[2][0] = Dz(1.,1.,h ) -2*Kex;
  tmp[3][0] = 0.;
  tmp[4][0] = 0.;
  tmp[5][0] = 0.;
  
  
  /* Parte de Intercambio  \/  */
  /*Coordenada (0,1)  y =1     */        row = &tmp[6*KX];        row[0][0]+= Kex; row[1][0]+= Kex; row[2][0]+= Kex;    
  /*Coordenada (0,-1) KY-1     */        row = &tmp[6*KX*(KY-1)]; row[0][0]+= Kex; row[1][0]+= Kex; row[2][0]+= Kex;    
  /*Coordenada (1,0)  x=1      */        row = &tmp[6];           row[0][0]+= Kex; row[1][0]+= Kex; row[2][0]+= Kex;    
  /*Coordenada (-1,0) KX-1     */        row = &tmp[6*(KX-1)];    row[0][0]+= Kex; row[1][0]+= Kex; row[2][0]+= Kex;    
  /* Parte de Intercambio  /\  */
  
  fprintf (stderr,".");
  fftw_execute(tensor);
  fprintf (stderr,".");
  fftw_destroy_plan (tensor); 
  
  /* Michele Voto :  Normalizar el tensor y no el campo por la linealidad de la FFT */
  for (int x=0;x<6*KX*KY;x++) 
  {
    nn[x][0] = (float) (tmp[x][0]*kvolume);
    nn[x][1] = (float) (tmp[x][1]*kvolume);
  }
  fftw_free(tmp);
  fftw_cleanup_threads();
  fprintf (stderr,".\n");
  return;
}


void savepng32 (int x,int y,void *ptr,char *filename)
{
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep * row_pointers;
  png_byte  *img;
  img = (png_byte*) ptr;

  FILE *fp = fopen(filename, "wb");
  if (!fp) return;
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) return ;
  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) return;
  if (setjmp(png_jmpbuf(png_ptr))) return;
  
  png_init_io(png_ptr, fp);

  if (setjmp(png_jmpbuf(png_ptr))) return ;

  png_set_IHDR(png_ptr,info_ptr,x,y,8,PNG_COLOR_TYPE_RGB_ALPHA,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_DEFAULT);
  png_set_compression_level(png_ptr,9); // Z_BEST_COMPRESSION
  
  png_write_info(png_ptr, info_ptr);
  
  if (setjmp(png_jmpbuf(png_ptr))) return;
  
  int row_bytes = png_get_rowbytes(png_ptr,info_ptr);
  row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * y);
  for (int row=0;row<y;row++) row_pointers[row] = (png_byte*) &(img[row*row_bytes]) ;
  
  png_write_image(png_ptr,row_pointers);
    
  if (setjmp(png_jmpbuf(png_ptr))) return;
  png_write_end(png_ptr, NULL);
  
  free(row_pointers);
  
  fclose(fp);
  return ;
}
