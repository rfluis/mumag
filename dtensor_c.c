#include <math.h>

// Replica reciente del codigo fortran 
// de Luis López Díaz

// Mayo de 2016 - Ricardo-Francisco Luis Martínez

#define pi     3.14159265358979323846264338327950288419716939937511
#define mfpi -12.56637061435917295385057353311801153678867759750042

double fnxx(double x,double y,double z,double dx,double dy,double dz);
double fnxy(double x,double y,double z,double dx,double dy,double dz);
static inline double f2(double x,double y,double z);
static inline double g2(double x,double y,double z);
static inline double f(double x,double y,double z);
static inline double g(double x,double y,double z);
static inline double fi (double x);
double Dz(double a,double b,double c);


double fnxx(double x,double y,double z,double dx,double dy,double dz)
{
  double answer;
  answer =  8.0 * f2(x,y,z) - 4.0 * ( f2(x,y,z+dz) + f2(x,y,z-dz) + f2(x,y+dy,z) + f2(x,y-dy,z) + f2(x+dx,y,z) + f2(x-dx,y,z) )
          + 2.0 * ( f2(x,y+dy,z+dz) + f2(x,y+dy,z-dz) + f2(x,y-dy,z+dz) + f2(x,y-dy,z-dz) + f2(x+dx,y+dy,z) + f2(x+dx,y-dy,z)
                  + f2(x-dx,y+dy,z) + f2(x-dx,y-dy,z) + f2(x+dx,y,z+dz) + f2(x+dx,y,z-dz) + f2(x-dx,y,z+dz) + f2(x-dx,y,z-dz) )
          -       ( f2(x+dx,y+dy,z+dz) + f2(x+dx,y+dy,z-dz) + f2(x+dx,y-dy,z+dz) + f2(x+dx,y-dy,z-dz) 
                  + f2(x-dx,y+dy,z+dz) + f2(x-dx,y+dy,z-dz) + f2(x-dx,y-dy,z+dz) + f2(x-dx,y-dy,z-dz) ) ;
      
      return answer / (dx*dy*dz*mfpi);
}
double fnxy(double x,double y,double z,double dx,double dy,double dz)
{
  double answer;
  answer = 8.0 * g2(x,y,z) - 4.0 * ( g2(x,y,z+dz) + g2(x,y,z-dz) + g2(x,y+dy,z) + g2(x,y-dy,z) + g2(x+dx,y,z) + g2(x-dx,y,z) )
         + 2.0 * ( g2(x,y+dy,z+dz) + g2(x,y+dy,z-dz) + g2(x,y-dy,z+dz) + g2(x,y-dy,z-dz) + g2(x+dx,y+dy,z) + g2(x+dx,y-dy,z)
                 + g2(x-dx,y+dy,z) + g2(x-dx,y-dy,z) + g2(x+dx,y,z+dz) + g2(x+dx,y,z-dz) + g2(x-dx,y,z+dz) + g2(x-dx,y,z-dz) )
        -        ( g2(x+dx,y+dy,z+dz) + g2(x+dx,y+dy,z-dz) + g2(x+dx,y-dy,z+dz) + g2(x+dx,y-dy,z-dz)
                 + g2(x-dx,y+dy,z+dz) + g2(x-dx,y+dy,z-dz) + g2(x-dx,y-dy,z+dz) + g2(x-dx,y-dy,z-dz) );

      return answer / (dx*dy*dz*mfpi);
}

static inline double f2(double x,double y,double z)
{
  return (f(x,y,z) - f(x,0.0,z) - f(x,y,0.0) + f(x,0.0,0.0));
}
static inline double g2(double x,double y,double z)
{
  return  g(x,y,z) - g(x,y,0.0);
}

static inline double f(double x,double y,double z)
{
  double x2,y2,z2,R,answer;
  x2 = x*x;
  y2 = y*y;   
  z2 = z*z;   
  R = sqrt(x2 + y2 + z2);   
  answer = 0.0;   
  
  if (z2>0.0)   
  {
    answer = (1.0/6.0) * (2*x2 - y2 - z2) * R;   
    if ((x2>0.0)&&(y2>0.0)) answer = answer - x * y * z * atan( y*z / (x*R));   
    if (y2>0.0)             answer = answer + 0.5 * y * (z2-x2) * fi(y / sqrt(x2+z2)); 
    if ((x2>0.0)||(y2>0.0)) answer = answer + 0.5 * z * (y2-x2) * fi(z / sqrt(x2+y2));  
  }   
  else    
  {
    answer = (1.0/6.0) * (2*x2 - y2) * R;   
    if ((x2>0.0)&&(y2>0.0)) answer = answer - 0.5 * y * x2 * fi(y/fabs(x));   
  }   
 return answer;
}
      
 static inline double g(double x,double y,double z)
 {
   double x2,y2,z2,R,answer;
   
   x2 = x*x;   
   y2 = y*y;
   z2 = z*z;  
   R = sqrt(x2 + y2 + z2);
   answer = - x * y * R / 3.0;
      
   if (z2>0.0)    
   {
    answer = answer - (z*z2/6.0) * atan( x * y / (z * R) ) + (y/6.0) * (3.0*z2 - y2) * fi(x/sqrt(y2+z2)) + (x/6.0) * (3.0*z2 - x2) * fi(y/sqrt(x2+z2));   
    if (y2>0.0)             answer = answer - (z*y2/2.0) * atan( x * z / (y * R) )  ;
    if (x2>0.0)             answer = answer - (z*x2/2.0) * atan( y * z / (x * R) )  ;
    if ((x2>0.0)&&(y2>0.0)) answer = answer + (x*y*z) * fi(z/sqrt(x2+y2)) ;
   }   
   else   
   {
    if (y2>0.0)             answer = answer - (y*y2/6.0) * fi(x/fabs(y));
    if (x2>0.0)             answer = answer - (x*x2/6.0) * fi(y/fabs(x));   
   }
   return answer;      
}

static inline double fi (double x)
{
  return ( log( x + sqrt( 1.0 + (x*x) )));
}

double Dz(double a,double b,double c)
{
      
  double a2,b2,c2,answer,sqrta2b2,sqrta2c2,sqrtb2c2,sqrta2b2c2;

  a2 = a*a;
  b2 = b*b;    
  c2 = c*c;    
      
  sqrta2b2   = sqrt(a2 + b2);
  sqrta2c2   = sqrt(a2 + c2);
  sqrtb2c2   = sqrt(b2 + c2);   
  sqrta2b2c2 = sqrt(a2 + b2 + c2);
      
  answer = -( ((b2-c2)/(2.0*b*c)) * log((sqrta2b2c2-a)/(sqrta2b2c2+a)) + ((a2-c2)/(2.0*a*c)) * log((sqrta2b2c2-b)/(sqrta2b2c2+b))
            + (b/(2.0*c))         * log((sqrta2b2+a)/(sqrta2b2-a)) 
            + (a/(2.0*c))         * log((sqrta2b2+b)/(sqrta2b2-b))
            + (c/(2.0*a))         * log((sqrtb2c2-b)/(sqrtb2c2+b))
            + (c/(2.0*b))         * log((sqrta2c2-a)/(sqrta2c2+a))
     + 2.0 * atan((a*b)/(c*sqrta2b2c2))
     + (a2*a+b2*b-2.0*c2*c) / (3.0*a*b*c)
     + (a2+b2-2.0*c2) / (3.0*a*b*c) * sqrta2b2c2
     + c/(a*b) * (sqrta2c2 + sqrtb2c2)
     - (pow(sqrta2b2,3.0)+pow(sqrtb2c2,3.0)+pow(sqrta2c2,3.0))/(3.0*a*b*c)) / pi;
  
  return answer;
      
}
