/*

Multiplier map
https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Mandelbrot_set#internal_coordinate_and_multiplier_map


algorithm: 
https://mathr.co.uk/blog/2013-04-01_interior_coordinates_in_the_mandelbrot_set.html
by Claude Heiland-Allen


code mostly based on the : 
 http://mathr.co.uk/blog/2014-11-02_practical_interior_distance_rendering.html
by 


 gcc -std=c99 -Wall -Wextra -pedantic -O3 -fopenmp b.c -lm
 
 
 
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double pi = 3.141592653589793;
const double infinity = 1.0 / 0.0;
const double colour_modulus = 5.7581917135421046e-2; // (1.0 + 1.0 / (phi * phi)) / 24.0;
const double ER2 = 2.0 * 2.0; // ER*ER
const double EPS2 = 1e-100 * 1e-100; // EPS*EPS





static inline double cabs2(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

static inline unsigned char *image_new(int width, int height) {
  return malloc(width * height * 3);
}

static inline void image_delete(unsigned char *image) {
  free(image);
}

static inline void image_save_ppm(unsigned char *image, int width, int height, const char *filename) {
  FILE *f = fopen(filename, "wb");
  if (f) {
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(image, width * height * 3, 1, f);
    fclose(f);
  } else {
    fprintf(stderr, "ERROR saving `%s'\n", filename);
  }
}

static inline void image_poke(unsigned char *image, int width, int i, int j, int r, int g, int b) {
  int k = (width * j + i) * 3;
  image[k++] = r;
  image[k++] = g;
  image[k  ] = b;
}

static inline void colour_hsv_to_rgb(double h, double s, double v, double *r, double *g, double *b) {
  double i, f, p, q, t;
  if (s == 0) { *r = *g = *b = v; } else {
    h = 6 * (h - floor(h));
    int ii = i = floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
    case 0: *r = v; *g = t; *b = p; break;
    case 1: *r = q; *g = v; *b = p; break;
    case 2: *r = p; *g = v; *b = t; break;
    case 3: *r = p; *g = q; *b = v; break;
    case 4: *r = t; *g = p; *b = v; break;
    default:*r = v; *g = p; *b = q; break;
    }
  }
}

static inline void colour_to_bytes(double r, double g, double b, int *r_out, int *g_out, int *b_out) {
  *r_out = fmin(fmax(255 * r, 0), 255);
  *g_out = fmin(fmax(255 * g, 0), 255);
  *b_out = fmin(fmax(255 * b, 0), 255);
}

static inline void colour_mandelbrot(unsigned char *image, int width, int i, int j, int period, double intRadius, double intAngle) {
  double r, g, b;
  
 /* if position > 1 then we have repetition of colors it maybe useful    */
      
  //if (intRadius > 1.0) intRadius -= 1.0; // normalize 
  
 
  
  
    
  colour_hsv_to_rgb(period * colour_modulus, 0.5, 1.0, &r, &g, &b); // changed b from tanh(distance )to 1.0 
  int ir, ig, ib;
  colour_to_bytes(r, g, b, &ir, &ig, &ib);
  
  //generate internal grid 
  if (intRadius <1.0){
    int rMax=10; /* number of color segments */
    double m=rMax* intRadius;
    int im=(int)m; // integer of m
 
    //if (intAngle > 1.0) intAngle -= 1.0; // normalize 
    int aMax=12; /* number of color segments */
    double k=aMax* intAngle;
    int ik =(int)k; // integer of m
    
     if ( im % 2 ) {ir -=40; }
     if ( ik % 2 ) {ig -=40; }
   }
   
   
   image_poke(image, width, i, j, ir, ig, ib);
   
}






/* newton function : N(z) = z - (fp(z)-z)/f'(z)) */

complex double N( complex double c, complex double zn , int pMax, double er2){

  
complex double z = zn;
complex double d = 1.0; /* d = first derivative with respect to z */
int p; 

for (p=0; p < pMax; p++){

   //printf ("p = %d ;  z = %f ; %f ;  d = %f ; %f \n", p, creal(z), cimag(z), creal(d), cimag(d)); 
   d = 2*z*d; /* first derivative with respect to z */
   z = z*z +c ; /* complex quadratic polynomial */
   //if (cabs(z) >er2) break;
    
}

 
 
 //printf (" next \n\n"); 
 //if ( cabs2(d) > 2) 
     z = zn - (z - zn)/(d - 1) ;
    
 

return z;
}



/* 
compute periodic point 
of complex quadratic polynomial
using Newton iteration = numerical method 

*/

complex double GivePeriodic(complex double c, complex double z0, int period, double eps2, double er2){

complex double z = z0;
complex double zPrev = z0; // prevoiuus value of z
int n ; // iteration
const int nMax = 64;


for (n=0; n<nMax; n++) {
     
    z = N( c, z, period, er2);
    if (cabs2(z - zPrev)< eps2) break;
    
    zPrev = z; }

return z;
}





complex double AproximateMultiplierMap(complex double c, int period, double eps2, double er2){
     
     complex double z;  // variable z 
     complex double zp ; // periodic point 
     complex double zcr = 0.0; // critical point
     complex double d = 1;
     complex double w;
     
     int p;
     
     zp =  GivePeriodic( c, zcr, period,  eps2, er2); // Find periodic point z0 such that f^p(z0,c)=z0 using Newton's method in one complex variable
     //zp = find(-1, 0, period, c); 
     //printf (" zp = %f ; %f p = %d \n", creal(zp), cimag(zp), period); 
     
     // Find w by evaluating first derivative with respect to z of f^p at z0 
     if ( cabs2(zp)<er2) {
     
     //printf (" zp = %f ; %f p = %d \n", creal(zp), cimag(zp), period); 
     z = zp;
     for (p=0; p < period; p++){
        d = 2*z*d; /* first derivative with respect to z */
        z = z*z +c ; /* complex quadratic polynomial */
     
          }
          
          
     
     
     }
     else d= 10000;

return d;
}




complex double MultiplierMap(complex double c, int period){

complex double w;
switch(period){
case 1: w = 1.0 - csqrt(1.0-4.0*c); 	break; // explicit
case 2: w = 4.0*c + 4; 			break; //explicit
default:w = AproximateMultiplierMap(c, period, EPS2, ER2); 	break; // 


}

return w;

}


// http://reference.wolfram.com/language/ref/MandelbrotSetBoettcher.html
complex double BoettcherMap(complex double c){


complex double c2 = c*c;
complex double c3 = c*c2;

return c + 0.5 - 1.0/(8.0*c) + 1.0/(16.0*c2) - 5.0/(128.0*c3) ; 
}



void ColorExterior(unsigned char *image, int width, int i, int j,double extRadius, double extAngle){

 double r, g, b;
  
 /* if position > 1 then we have repetition of colors it maybe useful    */
      
  //if (intRadius > 1.0) intRadius -= 1.0; // normalize 
  
 
  
  
    
  colour_hsv_to_rgb(colour_modulus, 0.4, 1.0, &r, &g, &b); // changed b from tanh(distance )to 1.0 
  int ir, ig, ib;
  colour_to_bytes(r, g, b, &ir, &ig, &ib);
  
  //generate external grid 
 // if (extRadius <1.0){
    int rMax=10; /* number of color segments */
    double m=rMax* extRadius;
    int im=(int)m; // integer of m
 
    //if (intAngle > 1.0) intAngle -= 1.0; // normalize 
    int aMax=20; /* number of color segments */
    double k=aMax* extAngle;
    int ik =(int)k; // integer of m
    
     if ( im % 2 ) {ir -=50; }
     if ( ik % 2 ) {ig -=50; }
   
   
   
   image_poke(image, width, i, j, ir, ig, ib);
}



/*
carg_t(z):=
 block(
 [t],
 t:carg(z)/(2*%pi),  // now in turns 
 if t<0 then t:t+1, // map from (-1/2,1/2] to [0, 1) 
 return(t)
)$

*/

double GiveTurn( double complex z){
double t;

  t =  carg(z);
  t /= 2*pi; // now in turns
  if (t<0.0) t += 1.0; // map from (-1/2,1/2] to [0, 1) 
  return (t);
}



static inline void render(unsigned char *image,  int pMax, int width, int height, complex double center, double radius) {

  double pixel_spacing = radius / (height / 2.0);

   #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < height; ++j) {
     for (int i = 0; i < width; ++i) {
      
      double x = i + 0.5 - width / 2.0;
      double y = height / 2.0 - j - 0.5;
      complex double c = center + pixel_spacing * (x + I * y);
      
      
      complex double w; // b(c):= (sqrt(1-4*c)+1)/4;
      double Radius;
      double Angle = 0.0;
      int period = 0;

      // find period and w 
      for (int p = 1; p < pMax; p++){
      
      if( cabs(c)<=2){
      w =  MultiplierMap(c,p);
      //if (p>2 ) printf (" w = %f ; %f p = %d \n", creal(w), cimag(w), p); 
      Radius = cabs(w);
      if ( Radius <= 1.0) {
           Angle = GiveTurn(w);
           period=p;
           //if (p>2 ) printf (" c = %f ; %f p = %d \n", creal(c), cimag(c), p); 
           break;}
         }   
      } 
      // colour 
       if (period==0) // exterior
         { w = BoettcherMap(c);
           Angle = GiveTurn(w); 
           Radius = cabs(w);
           ColorExterior(image, width, i, j,Radius, Angle);
           //image_poke(image, width, i, j, 0, 0, 0);
           }
       else colour_mandelbrot(image, width, i, j, period, Radius, Angle);
 
    }
  }
}

int main() {

  

  
 // int MaxIters = 100;
  int PeriodMax = 10;
  int width = 1000;
  int height = 1000;
  complex double center = -0.75 ;
  double radius = 1.5;
  const char *filename = "20.ppm";
 
  unsigned char *image = image_new(width, height);

  render(image,  PeriodMax, width, height, center, radius);
  image_save_ppm(image, width, height, filename);
  image_delete(image);
  
  
 

  return 0;
}
