/*
original algorithm by Bill Tavis posted on fractalforums.com
Re: smooth external angle of Mandelbrot set?
February 20, 2015, 05:12:34 PM
http://www.fractalforums.com/index.php?topic=19060.msg80726#msg80726

modified by Claude Heiland-Allen <claude@mathr.co.uk>
main program, bugfixes, collect bits when passing dwell bands
February 26, 2015
added initial bits from argument of first iterate
March 3, 2015

testing shows the original atan2() is only accurate to around 16 bits
bit collection when passing dwell bands is much more accurate

g++ -std=c++11 -Wall -Wextra -pedantic -O3 -o tavis tavis.cpp


#        cre   cim   iter  # pre(periodic) angle of nearby landing point
./tavis  0.251 0     1000  # .(0)
./tavis -2.001 0     1000  # .1(0)
./tavis -0.75  0.01  1000  # .(01)
./tavis  1e-16 1     1000  # .0(01)
./tavis -1.749 1e-10 1000  # .(011)
./tavis -1.01  0.251 1000  # .(01011001)
./tavis -0.99  0.251 1000  # .(01010110)

fails to terminate if command line arguments are invalid, for example:
  cre + cim i  not exterior to Mandelbrot set
  iter  too low for  cre + cim i  to escape
  
  
  
"  
... from every pixel, I drop a ball and let it roll downhill away from the set. 
Once it reaches a very large distance from the set, I use it's angle on that large circle. 
Every step of the ball rolling requires computing the gradient at its current location, 
which requires a full set of iterations. This is why it's so slow. 
However the advantage is that by decreasing the step size 
and using 4th order Runge-Kutta integration I was able to get it quite accurate:"

"  
  
  
  
*/

#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <vector>




typedef std::vector<bool> BITS;







BITS toBits(long double t) {
  BITS bs(64);
  t /= (8 * std::atan(1.0L));
  t -= std::floor(t);
  for (int i = 0; i < 64; ++i) {
    bs[i] = t >= 0.5;
    t *= 2;
    t -= floor(t);
  }
  return bs;
}






void printBits(const BITS &bs, size_t n) {
  std::cout << ".";
  for (size_t i = 0; i < bs.size(); ++i) {
    if (i == n) {
      std::cout << " ";
    }
    std::cout << (bs[i] ? "1" : "0");
  }
  std::cout << std::endl;
}




size_t bitsAgree(const BITS &as, const BITS &bs) {
  auto a = as.begin();
  auto b = bs.begin();
  size_t i = 0;
  bool a0 = false;
  bool b0 = false;
  while (a != as.end() && b != bs.end()) {
    a0 = *a++;
    b0 = *b++;
    i++;
    if (a0 != b0) { break; }
  }
  while (a != as.end() && b != bs.end()) {
    if (*a++ != b0 || *b++ != a0) { break; }
    i++;
  }
  return i;
}



/*
Here is a c++ function to get the gradient of a c-value based on the
formulas given by Linas Vepstas. 

imput : 
- c 
- maxiter
- 
output ( by reference ) :
- The gradient gz = gzx+gzy*i  
- the length 
- derivative Der 

to make it work with my other function below.

*/
void getGradient(long double cx, long double cy, long double &xDer,
                 long double &yDer, long double &gzx, long double &gzy,
                 long double &length, unsigned int maxIter,
                 long double &m, bool &b, long double &zarg) {
        // based on gradient equation found at:
        // http://linas.org/art-gallery/escape/ray.html
        // gradient = 2Df = f *zn*  Dzn / |zn|^2 * log |zn|


        
        
    long double zx = 0.0L;
    long double zy = 0.0L;
    xDer = 0;
    yDer = 0;
    long double xTemp;
    unsigned int n = 0;
    long double huge = std::numeric_limits<float>::max();
    while ((zx*zx + zy*zy) < huge && n < maxIter) {
            // compute new derivative:
        xTemp = 2.0L*(zx*xDer - zy*yDer) + 1.0L;
        yDer = 2.0L*(zx*yDer + zy*xDer);
        xDer = xTemp;
            // compute new z value:
        xTemp = zx*zx - zy*zy + cx;
        zy = 2.0L*zx*zy + cy;
        zx = xTemp;
        n++;
    }
    length = std::sqrt(zx*zx + zy*zy);
    zarg = std::atan2(zy, zx);
    b = zy < 0 || (zy == 0 && zx < 0);
    m = n + 1 - std::log2(std::log(length) / (0.5 * std::log(huge)));
    long double f = pow(2.0L,-m); // Douady-Hubbard potential
    
    
    long double Dx = zx*(xDer/2.0) - zy*(-yDer/2.0); // real part of zn * Dn
    long double Dy = zx*(-yDer/2.0) + zy*(xDer/2.0); // imag part of zn * Dn
    gzx = f*Dx / (length*length*std::log(length)); // x-gradient for cx,cy
    gzy = f*Dy / (length*length*std::log(length)); // y-gradient for cx,cy
}


// void getGradient(long double cx, long double cy, long double &xDer,
//                 long double &yDer, long double &gzx, long double &gzy,
//                 long double &length, unsigned int maxIter) {
//  long double m;
//  bool b;
//  long double zarg;
//  getGradient(cx, cy, xDer, yDer, gzx, gzy, length, maxIter, m, b, zarg);
// }











/*
And here is a function that uses the getGradeint function to perform
Runge-Kutta integration and find the external angle. It's assumed that
the point has already been checked for being outside the set.
Code:
*/

//BITS externalAngle(long double cx, long double cy, size_t maxIter, long double &e);

BITS externalAngle(long double cx, long double cy, size_t maxIter, long double &e) {
        // uses 4th order Runge-Kutta integration to find external angle
    long double ckx = cx;
    long double cky = cy;
    long double ckxTemp, ckyTemp;
    long double xDer,yDer,length,derLength,dist,gzx,gzy,gzxA,gzyA,gzxB,gzyB,gzxC,gzyC,gzxD,gzyD,gzlength;
    long double m = maxIter * 2;
    long double m_old = m;
    long double zarg = 0;
    bool first = true;
    BITS bs;
    bool b;
    while (m_old > 2) {
        // integrate with Runge-Kutta outwards along gradient:
        getGradient(ckx,cky,xDer,yDer,gzxA,gzyA,length,maxIter, m, b, zarg);
        
        
        if (first) {
          first = false;
          bs = toBits(zarg);
          reverse(bs.begin(), bs.end());
        }
        else
        if (std::ceil(m) < m_old) {
          bs.push_back(b);
          m_old = m;
        }
        derLength = std::sqrt(xDer*xDer + yDer*yDer);
        dist = 0.5L * std::log(length) * length / derLength;
        dist = std::min(std::max(16.0L, 0.01 * std::sqrt(ckx * ckx + cky * cky)),0.5*dist); // reduce step size
        gzlength = std::sqrt(gzxA*gzxA + gzyA*gzyA);
        ckxTemp = ckx + (0.5 * dist * gzxA/gzlength); // walk to midpoint
        ckyTemp = cky + (0.5 * dist * gzyA/gzlength); // using A
        getGradient(ckxTemp,ckyTemp,xDer,yDer,gzxB,gzyB,length,maxIter,m, b, zarg);
        gzlength = std::sqrt(gzxB*gzxB + gzyB*gzyB);
        ckxTemp = ckx + (0.5 * dist * gzxB/gzlength); // walk to midpoint
        ckyTemp = cky + (0.5 * dist * gzyB/gzlength); // using B
        getGradient(ckxTemp,ckyTemp,xDer,yDer,gzxC,gzyC,length,maxIter,m, b, zarg);
        gzlength = std::sqrt(gzxC*gzxC + gzyC*gzyC);
        ckxTemp = ckx + (dist * gzxC/gzlength); // walk full step
        ckyTemp = cky + (dist * gzyC/gzlength); // using C
        getGradient(ckxTemp,ckyTemp,xDer,yDer,gzxD,gzyD,length,maxIter,m, b, zarg);
        gzx = (1.0L/6.0L) * (gzxA + 2.0L*(gzxB + gzxC) + gzxD); // average the values together
        gzy = (1.0L/6.0L) * (gzyA + 2.0L*(gzyB + gzyC) + gzyD);
        gzlength = std::sqrt(gzx*gzx + gzy*gzy); // get length to normalize
        ckx += dist * gzx/gzlength; // normalize gradient and scale by dist
        cky += dist * gzy/gzlength;
    }
    // external angle
    e = (std::atan2(cky,ckx));
    e /= 2.0*M_PI; // now in turns
    if (e<0.0) e += 1.0; // map from (-1/2,1/2] to [0, 1) 
    //
    reverse(bs.begin(), bs.end());
    return bs;
}

/*
Even though Runge-Kutta method requires solving the gradient 4 times per
step, it winds up being much faster than Euler integration in practice
because it is so much more accurate. Notice that I used
Code:

dist = std::min(16.0L,0.5*dist);

for the step size. To get an image just as smooth with plain Euler
integration required reducing the step size to
Code:

dist = std::min(1.0L,0.05*dist);

which was way too slow.
*/



/*
Yes I am willing to but I don't have anything ready and it will take a little bit of work to extract something shareable... in the mean time I can share some pseudocode. The important function is one called gradientRay which is what calculates the point walking inward. This would be called from another function that is choosing angles to sample, and then plotting the orbits of the resulting c values. With simple Euler integration it looks like below. Larger values for the radius are more accurate in terms of walking down specific angles, but it takes longer. A similar thing can be said about reducing the step size, in terms of a trade off in accuracy vs time.

The equation to calculate the gradient can be found here  http://linas.org/art-gallery/escape/ray.html

Also more code can be found in this thread about external angles, which is doing the exact same thing but in reverse http://www.fractalforums.com/programming/smooth-external-angle-of-mandelbrot-set/

*/




// https://fractalforums.org/fractal-mathematics-and-new-theories/28/the-buddhas-jewel-revisited-finding-mandelbrot-orbits-shaped-like-julia-sets/2738/;topicseen
// Bill Tavis
function gradientRay(angle) (
     radius = 256;
     complex c; // declare a complex point
     c.real = radius * cos(angle);
     c.imag = radius * sin(angle);
     dist = radius; // declare distance variable
     stepsize = 0.5; // decrease the distance walked
     while (dist > 1.0e-19l) ( // walk in
          complex cG; // declare complex point for holding gradient
          getGradient(c,&cG,&dist) // function computes gradient and distance for the point c
          c -= stepsize * dist * (cG / cG.abs()); // take a step in the direction of the set
     )
     return c;
)




int main(int argc, char **argv) {


  if (argc > 3) {
    long double e = 0;
    auto bs2 = externalAngle(strtold(argv[1], 0), strtold(argv[2], 0), atoi(argv[3]), e);
    // test
    auto bs1 = toBits(e);
    auto n = bitsAgree(bs1, bs2);
    //printBits(bs1, n);
    
    printf("\ne = %.18Lf = ", e);
   // printBits(bs1, n);
    printBits(bs2, n);
    printf("\n");
    return 0;
  }
  
  printf(" not enough aruments : argc = %d and it should be argc > 3 \n", argc);
  
  return 1;
}
