
# How to compute external angle?



## Boettcher function
* [Boettcher equation in wikipedia](https://en.wikipedia.org/wiki/B%C3%B6ttcher%27s_equation)
* [Algorithm](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/boettcher#ArgPhi_-_External_angle_-_angular_component_of_complex_potential)


From the definitions:  
$`\Phi_M(c) \overset{def}{=} \Phi_c(c)`$  
$`arg_M(c)  = arg(\Phi_M(c)) `$  
$`arg_c(z) = arg(\Phi_c(z c)) `$  
so [external angle of point c on the parameter plane is equal to external angle of the point z=c on the dynamic plane](https://en.wikipedia.org/wiki/External_ray)  
$`arg_M(c) = arg_c(z=c) = arg(\Phi_c(z= c)) = arg(\Phi_M(c)) `$


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


## trace external ray outwards on the parameter plane

>>>
you need to trace a ray outwards, which means using different C values, and the bits come in reverse order, first the deepest bit from the iteration count of the start pixel, then move C outwards along the ray
(perhaps using the newton's method of mandel-exray.pdf in reverse), repeat until no more bits left.  you move C a fractional iteration count each time, and collect bits when crossing integer dwell boundaries
[Claude Heiland-Allen](http://mathr.co.uk/blog/)
>>>


##  Douady and Hubbard method for c near the real axis

>>>
Douady and Hubbard found a simple method for computing external angles for values of c outside of M and near the real axis. Call such an angle 2Pi*Ray, where 0 <= Ray < 1. The number Ray can be written as a binary decimal, i.e, as a sequence of zeroes and ones. 
To find it, consider the sequence  


    {Arg[c], Arg[c^2 +c], Arg[(c^2 + c)^2 + c], ...}.  
    

We replace Arg[z] by 
* 0 of 0 <= Arg[z] < Pi, 
* 1 otherwise.   

Here is some code for this.

c = -.75 +.0001*I; z = 0;

    Do[z = z^2 + c; Print[Abs[Floor[Arg[z]/Pi]]], {n, 1, 10}]

This produces the sequence {0, 1, 0, 1, 0, 1, 0, ...}, while c = -.75 - .0001*I produces {1, 0, 1, 0, 1, 0, 1, ...}. These are the binary expansions for 1/3 and 2/3 respectively. the point c0 = -.75 is the root of the period 2 bud. There are two rays leading inward to it, one coming from above and one from below. The two values of c we have chosen lie on or very near these two rays.
>>>
 [External angles in the Mandelbrot set: the work of Douady and Hubbard by Professor Douglas C. Ravenel](https://web.math.rochester.edu/people/faculty/doug/oldcourses/215s98/lecture10.html)








# See also
* [Computing external dynamic angle](https://gitlab.com/adammajewski/dynamic_external_angle)
* [parameter ray_in using Newton method and mpfr library in c ](https://gitlab.com/c_files/parameter_ray_in_newton_mpfr)
* [NonInteractive Parameter Ray_In MPFR](https://gitlab.com/adammajewski/NonInteractiveParameterRayInMPFR)
* [External angles in the Mandelbrot set: the work of Douady and Hubbard. by Professor Douglas C. Ravenel](https://web.math.rochester.edu/people/faculty/doug/oldcourses/215s98/lecture10.html)

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
* KaTeX to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html)

