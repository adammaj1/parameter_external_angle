
# How to compute external angle?

external angle 
* [external ray in wikipedia](https://en.wikipedia.org/wiki/External_ray)
* 


Methods:
* [series expansion formula of the Boettcher function](README.md#series-expansion-formula-of-the-boettcher-function)
* [trace external ray outwards on the parameter plane and collect bits ](README.md#trace-external-ray-outwards-on-the-parameter-plane) - the best method ?
* [Douady and Hubbard method](README.md#douady-and-hubbard-method-for-c-near-the-real-axis)
* [Stripe Average Coloring (or Method) = SAM or SAC](README.md#stripe-average-coloring-or-method-sam-or-sac) - good graphical result

## series expansion formula of the Boettcher function

### The Böttcher function
The Böttcher function maps the complement of the Mandelbrot set conformally to the complement of the closed unit disk.
* [Boettcher equation in wikipedia](https://en.wikipedia.org/wiki/B%C3%B6ttcher%27s_equation)
* [Algorithm](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/boettcher#ArgPhi_-_External_angle_-_angular_component_of_complex_potential)
* Böttcher function by Souichiro-Ikebe
  * [description](http://math-functions-1.watson.jp/sub1_spec_390.html#section060)
  * [How to process branch cut lines of Böttcher function](http://math-functions-1.watson.jp/sub4_math_020.html#section030)
  * [FractalRelated.m ](the Mathematica-Package file " FractalRelated.m " in "Code of special function graph" of "File created in Version 8 ") - package of the Mathematica code


From the definitions:  
$`\Phi_M(c) = \Phi_c(c)`$  
$`arg_M(c)  = arg(\Phi_M(c)) `$  
$`arg_c(z) = arg(\Phi_c(z)) `$  
so [external angle of point c on the parameter plane is equal to external angle of the point z=c on the dynamic plane](https://en.wikipedia.org/wiki/External_ray)  
$`arg_M(c) = arg_c(z=c) = arg(\Phi_c(z= c)) = arg(\Phi_M(c)) `$
[//]: # "\Phi_M(c) \overset{def}{=} \Phi_c(c)`"

where :
* $`arg_M`$ is exernal angle on the parameter plane
* $`arg_c`$ is the external angle on the dynamic plane
* $`arg`$ is the angle ( argument, phase) of complex number
* $`\Phi_c`$ is [the Boettcher map ](https://en.wikipedia.org/wiki/External_ray#Dynamical_plane_.3D_z-plane) on the dynamic plane
* $`\Phi_M`$ is [the Boettcher map ](https://en.wikipedia.org/wiki/External_ray#Dynamical_plane_.3D_z-plane) on the parameter plane


$` \Phi_c(z) = \lim_{n\to \infty} (f_c^n(z))^{2^{-n}} `$



Steps:
* on the parameter plane choose point c
* switch to the dynamic z-plane for c (not equal to zero)
  * take z= c 
  * check if it escapes = from exterior of filled Julia set
  * if yes: iterate z ( forward iteration  f_c(z) ) until it will be near infinity
  * remember number of final iteration when $`abs(z_n) > R_{limit} `$
  * switch plane to z-plane for c=0
* on the dynamic z-plane for c equal to zero
  * make n backward iterations 
  * angle of last z is aproximation of external angle 
  
  
  
Near infinity external angle of z is equal to angle of z so we can switch the plane


>" since the argument of the Boettcher function is not a main value but a value obtained by retaining multivalency, in fact, numerical calculation using this limit expression formula is very difficult. However, an infinite series easily obtained by transforming the limit expres>sion formula"  


Souichiro-Ikebe ( automatic translation)  



### series expansion formula for computing external angles

$`arg_M(c) = arg(c) + \sum_{n=1}^\infty \left( \frac{1}{2^n}*arg \left( \frac{f_c^n(c)}{f_c^n(c)-c}     \right ) \right )  `$


![whole set using palette colors](mturn.png)


![whole set using gray colors](mturng.png)


![zoom of wake 1/3](mturn3.png)


Code:

```c

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
double mturn(double complex c)
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
  
  
}//mturn
```



atan2(y/x) = Returns the principal value of the arc tangent of y/x, expressed in radians


How to compute arg(z/(z-c)) in a fast way?

simlify z/(z-c)

```c
z = x+y*I
arg(z) = atan2(y,x) 
z-c = u+v*I // b 
arg(z/(z-c)) = arg(z/b)
```

Here is Maxima CAS code:   
```
(%i2) z:x+y*%i;
(%o2)                                                                                                           %i y + x
(%i3) b:u+v*%i;
(%o3)                                                                                                           %i v + u
(%i4) creal(z/b);
                                                                                                                  %i y + x
(%o4)                                                                                                       creal(--------)
                                                                                                                  %i v + u
(%i5) realpart(z/b);
                                                                                                               v y + u x
(%o5)                                                                                                          ---------
                                                                                                                 2    2
                                                                                                                v  + u
(%i6) imagpart(z/b);
                                                                                                               u y - v x
(%o6)                                                                                                          ---------
                                                                                                                 2    2
                                                                                                                v  + u
(%i7) 

```

Denote :

$`z/b = rn/rd + I*in/id`$

notice that

$`atan2(in/id, rn/rd) = atan2(in, rn)`$

Now one can skip:
* 2 divisions : in/id and rn/rd
* 2 multiplications : $`v*v`$ and u*u
* 1 addition $`v^2+u^2`$





Files:
* [mturn.c](mturnc.c) - c file
* [mturn.png](mturn.png) - whole set using palette colors
* [mturng.png](mturng.png) - whole set using gray colors
* [mturn3.png](mturn3.png) - zoom of wake 1/3 using palette colors

Links:
* [description and cpp code by Wolf Jung](http://www.mndynamics.com/indexp.html#XR)
* [How to process branch cut lines of Böttcher function](http://math-functions-1.watson.jp/sub4_math_020.html#section030) by Souichiro-Ikebe
* [argphi](http://fraktal.republika.pl/cpp_argphi.html)



## trace external ray outwards on the parameter plane
>"mostly adopted a calculation method (of the external angle is) to trace the curve of the external ray" 

Souichiro-Ikebe ( automatic translation)  


>Testing shows the original atan2() is only accurate to around 16 bits, (so) bit collection when passing dwell bands is much more accurate.
```cpp
double externalAngle(...) {
...
	return (std::atan2(cy,cx));
}
```
>This gets you the angle in only double-precision, but using double precision floating point throughout it's possible to get the external angle in much higher precision 
>* the trick is to collect bits from the binary representation of the angle as you cross each dwell band 
>* whether the final iterate that escaped has a positive or negative imaginary part determines if the bit is 0 or 1 respectively, see [binary decomposition colouring](http://www.mrob.com/pub/muency/binarydecomposition.html)   

>You need to trace a ray outwards, which means using different C values, and the bits come in reverse order, first the deepest bit from the iteration count of the start pixel, then move C outwards along the ray
>(perhaps using [the newton's method of mandel-exray.pdf](http://www.math.titech.ac.jp/~kawahira/programs/mandel-exray.pdf) in reverse), repeat until no more bits left.  you move C a fractional iteration count each time, and collect bits when crossing integer dwell boundaries


Claude Heiland-Allen
* [blog](http://mathr.co.uk/blog/)
* [FF: smooth-external-angle-of-mandelbrot-set](http://www.fractalforums.com/programming/smooth-external-angle-of-mandelbrot-set/15/)  


See also:
* [binary-decomposition-and-external-angles](http://www.fractalforums.com/animations-showcase-(rate-my-short-animation)/binary-decomposition-and-external-angles/)
* [external-angle-contours](http://www.fractalforums.com/animations-showcase-(rate-my-short-animation)/external-angle-contours/)


Code:
* [tavis.cpp ](tavis.cpp) - compute external angle of point cx, cy

##  Douady and Hubbard method for c near the real axis

>
Douady and Hubbard found a simple method for computing external angles for values of c outside of M and near the real axis. Call such an angle 2Pi*Ray, where 0 <= Ray < 1.   
The number Ray can be written as a binary decimal, i.e, as a sequence of zeroes and ones. 
To find it, consider the sequence {Arg[c], Arg[c^2 +c], Arg[(c^2 + c)^2 + c], ...}.  


We replace Arg[z] by
* 0 if 0 <= Arg[z] < Pi, 
* 1 otherwise


Here is some Mathematica code for this.
```
    c = -.75 +.0001*I; 
    z = 0;
    Do[z = z^2 + c; Print[Abs[Floor[Arg[z]/Pi]]], {n, 1, 10}]
```   
This produces the sequence {0, 1, 0, 1, 0, 1, 0, ...} which is the binary expansion for 1/3  
For c = -.75 - .0001*I produces {1, 0, 1, 0, 1, 0, 1, ...} which is the binary expansions for 2/3.  
The point c0 = -.75 is the root of the period 2 bud. There are two rays leading inward to it, one coming from above and one from below. The two values of c we have chosen lie on or very near these two rays.
>  
   
[Douglas C. Ravenel](https://web.math.rochester.edu/people/faculty/doug/oldcourses/215s98/lecture10.html)


Files:
* [douady.c ](douady.c) - c file wich checks Douady-Hubbard method
* [morse.mac ](morse.mac) - batch file for Maxima cas which computes upper angles of external rays which land on the roots of the period doubling cascade on the real axis



## Stripe Average Coloring (or Method) = SAM or SAC


Links:
* [wikibooks](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/stripeAC)
* [gitlab](https://gitlab.com/adammajewski/mandelbrot_wiki_ACh)
* [wiki by Arnaoud Cheritat](https://www.math.univ-toulouse.fr/~cheritat/wiki-draw/index.php/Mandelbrot_set)
* [commons - whole set](https://commons.wikimedia.org/wiki/File:Mandelbrot_set_-_Stripe_Average_Coloring.png)
* [commons - wake 1/3](https://commons.wikimedia.org/wiki/File:Stripe_Average_Coloring_-_Mandelbrot_set_zoom_(_wake_1over3_).png)


```c
// the addend function
// input : complex number z
// output : double number t 
double Give_t(double complex z){

  return 0.5+0.5*sin(s*carg(z));

}

/*
  input :
  - complex number
  - intege
  output = average or other estimators, like distance or interior
 
*/
double Give_Arg(double complex C , int iMax)
{
  int i=0; // iteration 
   
   
  double complex Z= 0.0; // initial value for iteration Z0
  double A = 0.0; // A(n)
  double prevA = 0.0; // A(n-1)
  double R; // =radius = cabs(Z)
  double d; // smooth iteration count
  double complex dC = 0; // derivative
  double de; // = 2 * z * log(cabs(z)) / dc;
   
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    
      dC = 2 * Z * dC + 1; 
      Z=Z*Z+C; // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
      if (i>i_skip) A += Give_t(Z); // 
      
      R = cabs(Z);
      if(R > EscapeRadius) break; // exterior of M set
   
      prevA = A; // save value for interpolation
        
    } // for(i=0
   
   
  if (i == iMax) 
    A = -1.0; // interior 
  else { // exterior
    de = 2 * R * log(R) / cabs(dC);
    if (de < PixelWidth) A = FP_ZERO; //  boundary
    else {
      // computing interpolated average
      A /= (i - i_skip) ; // A(n)
      prevA /= (i - i_skip - 1) ; // A(n-1) 
      // smooth iteration count
      d = i + 1 + log(lnER/log(R))/M_LN2;
      d = d - (int)d; // only fractional part = interpolation coefficient
      // linear interpolation
      A = d*A + (1.0-d)*prevA;
     } 
  }
    
  return A;
}
```

and coloring: 
```c
  // compute 
  arg = Give_Arg( c, IterationMax);
     
  // 
  if (arg < 0.0)
     /*  interior of Mandelbrot set = inside_color =  */
      b = 0; 
    
  else // exterior and boundary 
     { 
     
     if (arg == FP_ZERO) b= 255; // boundary 
       else b = (unsigned char) (255 - 255*arg );
      
      
    };
    
    
    // gray gradient
      color[0]= b;  /* Red*/
      color[1]= b;  /* Green */ 
      color[2]= b;  /* Blue */
  return 0;
}
``` 

![whole set ](samm.png)



![zoom of wake 1/3](samm3.png)


Files:
* [samm.c](samm.c) - c file
* [samm.png](samm.png)
* [samm3.png](samm3.png)



# See also
* [Computing external dynamic angle](https://gitlab.com/adammajewski/dynamic_external_angle)
* [parameter ray_in using Newton method and mpfr library in c ](https://gitlab.com/c_files/parameter_ray_in_newton_mpfr)
* [NonInteractive Parameter Ray_In MPFR](https://gitlab.com/adammajewski/NonInteractiveParameterRayInMPFR)
* [External angles in the Mandelbrot set: the work of Douady and Hubbard. by Professor Douglas C. Ravenel](https://web.math.rochester.edu/people/faculty/doug/oldcourses/215s98/lecture10.html)
* [Plotting field lines during iteration by Chris Thomasson](http://www.fractalforums.com/new-theories-and-research/plotting-field-lines-during-iteration)


# git ( gitlab)

```
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/parameter_external_angle.git
git add .
git commit -m "Initial commit"
git push -u origin master
```


# technical note
GitLab uses:
* the Redcarpet Ruby library for [Markdown processing](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md)
* [KaTeX](https://khan.github.io/KaTeX/) to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html). [Here is used version](https://github.com/gitlabhq/gitlabhq/blob/a0715f079c143a362a7f6157db45020b8432003e/vendor/assets/javascripts/katex.js)

