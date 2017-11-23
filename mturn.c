
/* 



  "
  testing shows the original atan2() is only accurate to around 16 bits
  bit collection when passing dwell bands is much more accurate


  http://www.fractalforums.com/programming/smooth-external-angle-of-mandelbrot-set/15/
  
  double externalAngle(...) {
  ...
  return (std::atan2(cky,ckx));
}
  
  "This gets you the angle in only double-precision, but using double precision floating point throughout it's possible to get the external angle in much higher precision 
  - the trick is to collect bits from the binary representation of the angle as you cross each dwell band 
  - whether the final iterate that escaped has a positive or negative imaginary part determines if the bit is 0 or 1 respectively, 
  see binary decomposition colouring http://www.mrob.com/pub/muency/binarydecomposition.html .""
	Claude Heiland-Allen
  

  
  
   c console program:
   
   
   --------------------------------
   1. draws parameter plane for Fc(z)=z*z +c
   
   -------------------------------         
   2. technique of creating ppm file is  based on the code of Claudio Rocchini
   http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
   create 24 bit color graphic file ,  portable pixmap file = PPM 
   see http://en.wikipedia.org/wiki/Portable_pixmap
   to see the file use external application ( graphic viewer)
   
   
   complex point c -> virtual 2D array -> memory 1D array -> ppm file on the disc -> png file 
   
   C -> pixel (iX,iY)  -> index k  -> 24bit color 
   
   -----
   https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
   complex numbers are built in type 
 
   --------------
   formated with emacs
   -------------
   to compile : 

 
 
   gcc mturn_ppm.c -lm -Wall 
 
 
   ./a.out
   
   
   to convert to png using ImageMagic

   convert phase.ppm phase.png  



  ----------------------
	



 
 
*/
#include <stdio.h>
#include <stdlib.h>		// malloc
#include <math.h>
#include <complex.h> // https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
 
 

 


 
/* screen ( integer) coordinate */

const int iWidth  = 1000; 
const int iHeight = 1001; // +1 for the main antenna 


/* world ( double) coordinate = parameter plane*/
// double complex C =  Cx + Cy*I ;
double CxMin;
double CxMax;
double CyMin;
double CyMax;

/* */
double PixelWidth; //=(ZxMax-ZxMin)/iWidth;
double PixelHeight; // =(ZyMax-ZyMin)/iHeight;

double thickness = 0.50;// multiplier 
double MinBoundaryWidth; //   


// ------------ colors ----------------

/* color component ( R or G or B) is coded from 0 to 255 */
/* it is 24 bit color RGB file */
int ColorBytes = 3; // 3*8 = 24 bit color 

#define iBlack  0 
#define iWhite  255 
#define iRed 	1 
#define iGreen  2   
#define iBlue  	3  
#define iCyan  	4   

int repetition= 150; // of colors

/* -------  iterations  */
#define IterationMax 40 

/* bail-out value , radius of circle ;  */
const double EscapeRadius=2.0;
double m = 2.0; //  multiplier


double TwoPi=2.0*M_PI;

 


//  memory 1D array  = for the image processing
unsigned char *data;       
size_t MemmorySize;   


// ---------------- color functions -------------------------------
       
void GiveGrayColor(double position , int k, unsigned char c[]) 
{
  unsigned char b = 255*position;
  c[k]   = b;
  c[k+1] = b;
  c[k+2] = b;

}


void ColorPixel(int iColor, int k, unsigned char c[])
{
  switch (iColor)
    {
    case iBlack : 	c[k]   = 0; 	c[k+1] = 0;	c[k+2] = 0; break;
    case iRed:    	c[k]   = 255; 	c[k+1] = 0;	c[k+2] = 0;   break;
    case iGreen : 	c[k]   = 0; 	c[k+1] = 255;	c[k+2] = 0; break;
    case iBlue : 	c[k]   = 0;	c[k+1] = 0;	c[k+2] = 255; break;
    case iCyan : 	c[k]   = 0; 	c[k+1] = 255;	c[k+2] = 255; break;
    case iWhite : 	c[k]   = 255; 	c[k+1] = 255;	c[k+2] = 255; break;
    }
}

// palette ; when position > 1 then there is repetition of colors
void GiveMandelColor(double position, int k, unsigned char c[]) 
{
  double b = position - floor(position);// only fractional part = repetition of colors when position >1
  
  if (b< 0.25){ColorPixel(iRed, k, c); return ;} 
  if (b< 0.50){ColorPixel(iGreen, k, c); return ;} 
  if (b< 0.75){ColorPixel(iBlue, k, c); return ;} 
  //if (b< 1.00)
  ColorPixel(iCyan, k, c); return ; 

}




       
/* 
   gives position ( index) in 1D virtual array  of 2D point (iX,iY) from ; uses also global variable iWidth 
   without bounds check !!
*/
int f(int ix, int iy)
{ return ColorBytes*(ix + iy*iWidth); }
        
        
        
double complex give_c(int iX, int iY){
  double Cx, Cy;
  Cy = CyMax - iY*PixelHeight; // inverse y axis
  
  Cx = CxMin + iX*PixelWidth;
   
  return Cx + Cy*I;
 
 
}
 
 

int IsInside(double complex c)
{
  double complex z = 0.0;
  int i;
  int iMax = 400;
  
  for (i=0; i<iMax; i++){
  	if (cabs(z)>EscapeRadius) break;
  	z= z*z+c;
  }
  
  if (i==iMax) return 1; // is Inside
     
  return 0 ; // outside
}



/*
  input :
  - complex number
  - intege
  output =  estimator dn
   
*/
double Give_D(double complex C , int iMax)
{
  int i=0; // iteration 
   
   
  double complex Z= 0.0; // initial value for iteration Z0
  double R; // =radius = cabs(Z)
  double D; // output of the function
  double complex dC = 0; // derivative with respect to c 
  // dn in A Cheritat notation
  double de; // = 2 * z * log(cabs(z)) / dc =  estimated distance to the Mandelbrot set
   
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    
      dC = 2 * Z * dC + 1; 
      Z=Z*Z+C; // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
            
      R = cabs(Z);
      if(R > EscapeRadius) break; // exterior of M set
   
      
    } // for(i=0
   
   
  if (i == iMax) D = -1.0; // interior 
    else { // exterior
      de = 2.0 * R * log(R) / cabs(dC) ; // 
    
      if (de < MinBoundaryWidth) D = FP_ZERO; //  boundary
         else  D = 1.0; // exterior
    }
    
  return D;  
}





//------------  function mturn520   ----------------------------------------------

/*
 Fc(z) = z*z + c
 z= x+y*i
 c= a+b*i
 This function computes external argument of point C in turns
 for mandelbrot set for Fc(z)= z*z + c
 external argument = Arg(Phi(c))
 1 [turn] = 360 [degrees] = 2* M_PI [radians]
 this function is based on function mturn from mbrot.cpp 
 from old (probably 5.2 of May 17, 2008) version of program mandel by Wolf Jung
 http://www.iram.rwth-aachen.de/~jung/indexp.html 

Already checked that escaping. Requires z = c instead of z = 0
 http://fraktal.republika.pl/cpp_argphi.html
 algorithm: http://www.mndynamics.com/indexp.html#XR

this formula will be valid not only for large |z| 

*/
 
double mturn520(double complex c)
{ 
  int j; 
  int jMax = 100;
  
  double s = 1.0; // = 1/(2^n)
  double dr = 1.0/2.0; 
  double theta; 
  // z-c 
  double u; // creal(z-c)
  double v; // cimag(z-c)
  
  double r; //   r here means r^2 = cabs(z)* cabs(z)
  
  // z = x+y*I 
  double x; 
  double y;
  
  // c = cx + cy*I        
  double cx = creal(c);
  double cy = cimag(c);        
  
  
  
  
  // Requires z = c instead of z = 0
  x = cx;
  y = cy; 
  
  // theta = arg(z)
  theta = atan2(y, x);
  
  // compute the sum 
  for (j = 1; j < jMax; j++)
  { 
    s *= dr;
    // z=Fc(z)
    double temp = x*x - y*y + cx;
    y = 2*x*y + cy;
    x = temp; 
    
    // r here means r^2 = cabs(z)* cabs(z) 
    r = x*x  + y*y; // 
    if (r < .0001) return -7; // ?
    
    
    // z-c
    u = x - cx; 
    v = y - cy;
    
    
    
    // atan2 here is computing  arg(z/(z-c))
    // s*atan2 means : theta/(2^n)
    // theta += is summation
    theta += s * atan2(u*y - v*x, u*x + v*y);
    //
    if (r/s > 1e25) break; // prevent -nan
  }
  
  
  //if (r < 1000) return -6; // ? lazy escaping ???, 
  //
  theta *= (.5/M_PI); // convertion to turns 
  //
  theta -= floor(theta); // modulo 1 
  
  return theta;
  
  
}//mturn520



// compute and save to the data array color of 1 pixel
int ComputeAndSavePixelColor(int iX, int iY){
 
  
  
  
  double complex C; // 
  int k; // index of the 1D array
  
  double t; // argument of complex numer in turns
  
     
  
  
  // index of 1D memory array
  k = f(iX, iY);   
  
  C = give_c(iX, iY);
  
  if (IsInside(C)) 
  	ColorPixel(iWhite, k, data);
  	else{ // exterior or boundary
  		// DEM/M
		double d = Give_D(C, 1000);
		if (d == FP_ZERO) 
			ColorPixel(iBlack, k, data); // boundary
			else{// exterior = external angle ( argument)	
				t =  mturn520(C);
  				//GiveGrayColor(t,k,data);
  				GiveMandelColor(t*repetition, k, data);
  				}  	
  		
  	}
  
  
  
    
 
   
  return 0;
}
 
 
 
int setup(){

  // whole set
  CxMin= -2.2;
  CxMax=  0.8;
  CyMin= -1.5;
  CyMax=  1.5;


  // bulb 3 zoom    
  // c = -0.095000000000000  +0.990000000000000 i    period = 0
  // c = -0.284736659610103  +1.001858541225631 i    period = 0
  /*
  CxMin= -0.095-0.189;
  CxMax= -0.095+0.189;
  CyMin=  1.0-0.189;
  CyMax=  1.0+0.189;
  */
  
  //
  PixelWidth=(CxMax-CxMin)/iWidth;
  PixelHeight=(CyMax-CyMin)/iHeight;
  //
   MinBoundaryWidth = thickness*PixelWidth;  
  //
  MemmorySize = iWidth * iHeight * ColorBytes * sizeof (unsigned char);	// https://stackoverflow.com/questions/492384/how-to-find-the-sizeof-a-pointer-pointing-to-an-array
        
  /* create dynamic 1D arrays for colors   */
  data = malloc (MemmorySize);
  if (data == NULL )
    { fprintf (stderr, " Error: Could not allocate memory !!! \n"); return 1; }

  printf (" No errors. End of setup \n");
  return 0;

}
 
 
 
 
 
// save dynamic "A" array to pgm file 
int SaveArray_2_PPM_file (unsigned char A[])
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char *filename = "mturn514_whole_mandel_150.ppm";
  char *comment = "# parameter plane";		/* comment should start with # */

  /* save image to the pgm file  */
  fp = fopen (filename, "wb");	/*create new file,give it a name and open it in binary mode  */
  fprintf (fp, "P6\n %s\n %u %u\n %u\n", comment, iWidth, iHeight, MaxColorComponentValue);	/*write header to the file */
  fwrite (A, MemmorySize, 1, fp);	/*write A array to the file in one step */
  printf ("File %s saved. \n", filename);
  fclose (fp);

  return 0;
}


 
 
void CreateImage(){
  int iX,iY; // screen = integer coordinate in pixels       

  // fill the array = render image = scanline 2D  of virtual 2D array 
  for(iY=0;iY<iHeight;iY++)
    for(iX=0;iX<iWidth;iX++)
      ComputeAndSavePixelColor(iX, iY); 
      	
      	
  SaveArray_2_PPM_file (data);     	  
} 
 
 
 
void info(){

  printf(" Parameter plane ( c plane) Mandelbrot set for complex quadratic polynomial fc(z) = z^2 + c\n ");
  printf(" Rectangle part of 2D parmeter plane: corners: \n CxMin = %f;   CxMax = %f;  CyMin = %f; CyMax = %f \n ", CxMin, CxMax, CyMin, CyMax);
  printf(" center and radius: \n CenterX = %f;   CenterY = %f;  radius = %f\n ", (CxMax+CxMin)/2.0, (CyMax+CyMin)/2.0, fabs(CyMax-CyMin)/2.0);
  printf(" Mag = zoom = %f\n ",  2.0/fabs(CyMax-CyMin));
  printf(" PixelWidth = %f and PixelHeight =%f\n", PixelWidth, PixelHeight);
  printf(" Escape Radius = %f\n ", EscapeRadius);
  printf(" Iteration Max = %d\n ", IterationMax);
  // printf(" M_PI = %.16f ; TwoPi = %.16f \n", M_PI, TwoPi);


} 
 
 
 
void close(){
 
  info(); 
  free(data); 
}
 
 
 
int main()
{
 
  setup();      
  CreateImage();     
  close();
  
        
  return 0;
}
