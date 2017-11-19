/*
gcc douady.c -Wall -lm

./a.out
 
  
*/




#include <stdio.h>
#include <math.h>
#include <complex.h> // 


/* 
   ========================================
   Bit Operations as macros
   ======================================== 
   http://www.mathcs.emory.edu/~cheung/Courses/255/Syllabus/1-C-intro/bit-array.html
   
   
  Map the bit array onto the array of integer as follows: 

 
 Shun Yan Cheung,
 A is array of int 
 int A[LENGTH];
   
*/
#define SetBit(A,k)     ( A[(k/32)] |= (1 << (k%32)) ) // set bit to 1 
#define ClearBit(A,k)   ( A[(k/32)] &= ~(1 << (k%32)) ) // set bit to 0 
#define TestBit(A,k)    ( A[(k/32)] & (1 << (k%32)) ) // if bit is 1 then true else false 

#define LENGTH 2  // length of the array with 32-bit cells 
int iMax = LENGTH*sizeof(int)*8; // length of the array in bits = maximal iteration !!!

double TwoPi=2.0*M_PI;


void ClearArray(int A[]){
  int l; // number of the array cell 
  int lMax = LENGTH;
  
  for ( l = 0; l < lMax; l++ )
      A[l] = 0;                    // Clear the bit array

 }
 
 
void PrintBinaryFraction(int A[]){

  int i; // bit number = iteration number 
  //      
   printf(";\t 0."); // arg(phi(c))     
   for ( i = 0; i < iMax; i++ )
      if ( TestBit(A, i ) )
         printf("1");
         else printf("0"); 
   printf("\n");     

 }
 
 
 void PrintComplex(double complex c){
 
 printf("c = (%.16f ; %.16f )", creal(c), cimag(c));
 }
 
 
 
 double GiveTurn(double complex z)
{
  double argument;
 
  argument = carg(z); //   argument in radians from -pi to pi
  if (argument<0) argument=argument + TwoPi; //   argument in radians from 0 to 2*pi
  return argument/TwoPi ; // argument in turns from 0.0 to 1.0
}

 
// https://web.math.rochester.edu/people/faculty/doug/oldcourses/215s98/lecture10.html

  
 int iterate(double complex C )
  {
   int i; // bit number = iteration number 
   // array of int as a array of bits !!
   int A[LENGTH]={0}; // set all bits to 0
   double complex Z = C; // initial value for iteration Z0
   
   
   
   for(i=0;i<iMax;i++)
    {
      if (GiveTurn(Z)>0.5)   SetBit( A, i ); // 1 
      Z=Z*Z+C; //       
    }
   PrintComplex(C); 
   PrintBinaryFraction(A);
    
   return 0; 
 }

/*

./tavis  0.251 0     1000  # .(0)  = decimal 0.0
./tavis -2.001 0     1000  # .1(0) = decimal 0.5
./tavis -0.75  0.01  1000  # .(01) = decimal 1/3 = 0.(3)
        -.75  -0.0001         .(10)  = {1, 0, 1, 0, 1, 0, 1, ...} = decimal 1/3 
./tavis  1e-16 1     1000  # .0(01) = decimal 1/6 = 0.1(6) = 0.1666666666666666... 
./tavis -1.749 1e-10 1000  # .(011)
./tavis -1.01  0.251 1000  # .(01011001)
./tavis -0.99  0.251 1000  # .(01010110)




iMax = 64 
it should c be from exterior of Mandelbrot set and near real axis 
http://fraktal.republika.pl/mset_external_ray_f1_2.html

c near cusp of main cardioid: rays 0 and 1
c = (0.2510000000000000 ; 0.0100000000000000 );	 0.0000000000000000000000000000000000000000000000000000000000000000
c = (0.2510000000000000 ; -0.0100000000000000 );	 0.1111111111111111111111111111111111111111111111111111111111111111

c near the root of period 2 component: rays 1/3= 0.(01) and 2/3= 0.(10)
c = (-0.7500000000000000 ; 0.0010000000000000 );	 0.0101010101010101010101010101010101010101010101010101010101010101
c = (-0.7500000000000000 ; -0.0010000000000000 );	 0.1010101010101010101010101010101010101010101010101010101010101010

c near the root of period 4 component: rays (6/15 = 2/5 = 0.(0110) , 9/15 = 0.(1001)) 
c = (-1.2500000000000000 ; 0.1000000000000000 );	 0.0100011001000101100110110111011000000000000000000000000000000000
c = (-1.2500000000000000 ; -0.1000000000000000 );	 0.1011100110111010011001001000100000000000000000000000000000000000

c near the root of period 8 component: rays (6/15 = 2/5 = 0.(0110) , 9/15 = 0.(1001)) 
c = (-1.3680989394000000 ; 0.0010000000000000 );	 0.0110000101100001011000010110000101100001011010010110100101101001
c = (-1.3680989394000000 ; -0.0010000000000000 );	 0.1001111010011110100111101001111010011110100101101001011010010110
c = (0.0000000000000001 ; 1.0000000000000000 );	 0.0010101010101010101010101010101010101010101100101110101000000000
c = (-1.7490000000000001 ; 0.0000000001000000 );	 0.0110110110110110110110110110110110110110110100100100100100100100
c = (-1.0100000000000000 ; 0.2510000000000000 );	 0.0101010101010101010101010101010101010101010101010101010101010101
c = (-0.9900000000000000 ; 0.2510000000000000 );	 0.0101010101010101010101010101010101010101010101010101010101010101

c near the end of the main antenna: ray 1/2
c = (-2.0009999999999999 ; 0.0001000000000000 );	 0.0111111111110010000000000000000000000000000000000000000000000000
c = (-2.0009999999999999 ; -0.0001000000000000 );	 0.1000000000001100000000000000000000000000000000000000000000000000



*/
int main(  )
{   
   
   printf("iMax = %d \n", iMax);
   printf("it should c be from exterior of Mandelbrot set and near real axis \n");
   printf("http://fraktal.republika.pl/mset_external_ray_f1_2.html\n");
   
   printf("\nroots on the real axis\n");
   printf("\nc near cusp of main cardioid: rays 0 and 1\n");
   iterate(0.251+0.01*I);
   iterate(0.251-0.01*I);
   
   printf("\nc near the root of period 2 component: rays 1/3= 0.(01) and 2/3= 0.(10)\n");
   iterate(-0.75 +0.001*I);
   iterate(-0.75 -0.001*I);
   printf("\nc near the root of period 4 component: rays (6/15 = 2/5 = 0.(0110) , 9/15 = 0.(1001)) \n");
   iterate(-1.25 +0.1*I);
   iterate(-1.25 -0.1*I);
   printf("\nc near the root of period 8 component: rays (6/15 = 2/5 = 0.(0110) , 9/15 = 0.(1001)) \n");
   iterate( -1.3680989394 +0.001*I);
   iterate( -1.3680989394 -0.001*I);
   
   
   iterate(1e-16+1.0*I);
   iterate(-1.749+1e-10*I);
   iterate(-1.01+0.251*I);
   iterate(-0.99+0.251*I);

   printf("\nc near the end of the main antenna: ray 1/2\n");
   iterate(-2.001+0.0001*I);
   iterate(-2.001-0.0001*I);

   
           
         
  return 0;       
         
}
