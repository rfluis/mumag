#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
/* Para agilizar las llamadas , inverse, jump y n serán variables globales
 * y son parametros de doFFT. En los lotes de FFT multidimensional son siempre
 * las mismas, cambiando solo el puntero de acceso
 */

int n,jump,inverse;

void fft3ds(int x,int y,int z,int inv,int s,double *GRe, double *GIm)
{
    /* Paso 1: Salto = 1 . Tamaño = Nz . Avance = Nz·y+Nz·Ny·x
     * Paso 2: Salto = Nz . Tamaño = Ny . Avance = z + Nz·Ny·x
     * Paso 3: Salto = Nz·Ny . Tamaño = Nx . Avance = z + Nz·y
     * 
     * for(i=0;i<x;i++) for(j=0;j<y;j++) for(k=0;k<z;k++)
     */
    inverse=inv;
    int i,j,k;
    
    jump=s*1; n=z;
    for(i=0;i<x;i++) for(j=0;j<y;j++) doFFT(GRe+s*(z*j+z*y*i),GIm+s*(z*j+z*y*i));
    jump=s*z; n=y;
    for(i=0;i<x;i++) for(k=0;k<z;k++) doFFT(GRe+s*(k+z*y*i),GIm+s*(k+z*y*i));
    jump=s*y*z; n=x;
    for(j=0;j<y;j++) for(k=0;k<z;k++) doFFT(GRe+s*(k+z*j),GIm+s*(k+z*j));
}

void fft3d(int x,int y,int z,double *GRe, double *GIm,int inv) 
{
    /* Paso 1: Salto = 1 . Tamaño = Nz . Avance = Nz·y+Nz·Ny·x
     * Paso 2: Salto = Nz . Tamaño = Ny . Avance = z + Nz·Ny·x
     * Paso 3: Salto = Nz·Ny . Tamaño = Nx . Avance = z + Nz·y
     * 
     * for(i=0;i<x;i++) for(j=0;j<y;j++) for(k=0;k<z;k++)
     */
    inverse=inv;
    int i,j,k;
    
    jump=1; n=z;
    for(i=0;i<x;i++) for(j=0;j<y;j++) doFFT(GRe+z*j+z*y*i,GIm+z*j+z*y*i);
    jump=z; n=y;
    for(i=0;i<x;i++) for(k=0;k<z;k++) doFFT(GRe+k+z*y*i,GIm+k+z*y*i);
    jump=y*z; n=x;
    for(j=0;j<y;j++) for(k=0;k<z;k++) doFFT(GRe+k+z*j,GIm+k+z*j);
}

void fft2ds(int x,int y,int inv,int s,double *GRe, double *GIm)
{
    int index;
    inverse=inv;
    // Paso 1
    // x transformadas de tamaño y
    jump=s;
    n=y;
    for(index=0;index<x;index++) doFFT(GRe+s*index*y,GIm+s*index*y);
    jump=s*y;
    n=x;
    for(index=0;index<y;index++) doFFT(GRe+s*index,GIm+s*index);
return;    
}

void fft2d(int x,int y,int inv,double *GRe, double *GIm)
{
    int index;
    inverse=inv;
    // Paso 1
    // x transformadas de tamaño y
    jump=1;
    n=y;
    for(index=0;index<x;index++) doFFT(GRe+index*y,GIm+index*y);
    jump=y;
    n=x;
    for(index=0;index<y;index++) doFFT(GRe+index,GIm+index);
return;    
}

void doFFT(double *GRe, double *GIm)
{
  //Calculate m=log_2(n)
  int m = ilog(n);
  int shift = 32 - m;
  int j,i,i1,l1,l2,l;
  double u1,u2,t1,t2,ca,sa;
  register double z;
  
  for(i=1;i<n;i++)
  {
      j=reverse_integer(i<<shift);
      if (j>i)
      {
          z     =GRe[j*jump];
          GRe[j*jump]=GRe[i*jump];
          GRe[i*jump]=z;
          
          z     =GIm[j*jump];
          GIm[j*jump]=GIm[i*jump];
          GIm[i*jump]=z;
      }
      
  }
  
    //Calculate the FFT
  ca = -1.0;
  sa = 0.0;
  l1 = 1;
  l2 = 1;
  
  for(l = 0; l < m; l++)
  {
    l1 = l2;
    l2 *= 2;
    u1 = 1.0;
    u2 = 0.0;
    for(j = 0; j < l1; j++)
    {
      for(i = j; i < n; i += l2)
      {
        i1 = i + l1;
        t1 = u1 * GRe[i1*jump] - u2 * GIm[i1*jump];
        t2 = u1 * GIm[i1*jump] + u2 * GRe[i1*jump];
        GRe[i1*jump] = GRe[i*jump] - t1;
        GIm[i1*jump] = GIm[i*jump] - t2;
        GRe[i*jump] += t1;
        GIm[i*jump] += t2;
      }
      z  = u1 * ca - u2 * sa;
      u2 = u1 * sa + u2 * ca;
      u1 = z;
    }
    if(inverse) sa = -sqrt((1.0 - ca) * 0.5); 
    else sa = sqrt((1.0 - ca) * (0.5));
    ca = sqrt((1.0 + ca) * 0.5);
  }
}

void fft1d(int n,int inverse, double *GRe, double *GIm)
{
  //Calculate m=log_2(n)
  int m = ilog(n);
  int shift = 32 - m;
  int j,i,i1,l1,l2,l;
  double u1,u2,t1,t2,ca,sa;
  register double z;
  
  for(i=1;i<n;i++)
  {
      j=reverse_integer(i<<shift);
      if (j>i)
      {
          z     =GRe[j];
          GRe[j]=GRe[i];
          GRe[i]=z;
          
          z     =GIm[j];
          GIm[j]=GIm[i];
          GIm[i]=z;
      }
      
  }
  
    //Calculate the FFT
  ca = -1.0;
  sa = 0.0;
  l1 = 1;
  l2 = 1;
  
  for(l = 0; l < m; l++)
  {
    l1 = l2;
    l2 *= 2;
    u1 = 1.0;
    u2 = 0.0;
    for(j = 0; j < l1; j++)
    {
      for(i = j; i < n; i += l2)
      {
        i1 = i + l1;
        t1 = u1 * GRe[i1] - u2 * GIm[i1];
        t2 = u1 * GIm[i1] + u2 * GRe[i1];
        GRe[i1] = GRe[i] - t1;
        GIm[i1] = GIm[i] - t2;
        GRe[i] += t1;
        GIm[i] += t2;
      }
      z  = u1 * ca - u2 * sa;
      u2 = u1 * sa + u2 * ca;
      u1 = z;
    }
    if(inverse) sa = -sqrt((1.0 - ca) * 0.5); 
    else sa = sqrt((1.0 - ca) * (0.5));
    ca = sqrt((1.0 + ca) * 0.5);
  }
}

unsigned int reverse_integer (unsigned int v)
{
  unsigned int c;
  const unsigned int table[256] ={
  0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
  0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
  0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
  0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
  0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
  0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
  0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
  0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
  0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
  0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
  0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
  0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
  0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
  0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
  0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
  0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};
  return (table[v & 0xff] << 24) | (table[(v >> 8) & 0xff] << 16) | (table[(v >> 16) & 0xff] << 8) | (table[(v >> 24) & 0xff]);
}

unsigned int ilog (unsigned int x)
{
  unsigned int ans=0;
  while (x=x>>1) ans++;
  return ans;
}


