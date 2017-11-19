
# How to compute external angle?



## Boettcher function
* [Boettcher equation in wikipedia](https://en.wikipedia.org/wiki/B%C3%B6ttcher%27s_equation)
* [Algorithm](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/boettcher#ArgPhi_-_External_angle_-_angular_component_of_complex_potential)


[External angle of point c on the parameter plane is equal to external angle of the point z=c on the dynamic plane](https://en.wikipedia.org/wiki/External_ray)
$`arg_M(c) = arg_c(z=c) = arg(\Phi_c(z= c)) = arg(\Phi_M(c)) `$


where :
* $`arg_M`$ is exernal angle on the parameter plane
* $`arg_c`$ is the external angle on the dynamic plane
* $`arg`$ is the angle ( argument, phase) of complex number
* $`\Phi_c`$ is [the Boettcher map ](https://en.wikipedia.org/wiki/External_ray#Dynamical_plane_.3D_z-plane)


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
  
  
  
Near infinity angle of z is equal to external angle of z so we can switch the plane


## trace external ray outwards on the parameter plane

>>>
you need to trace a ray outwards, which means using different C values,
and the bits come in reverse order, first the deepest bit from the
iteration count of the start pixel, then move C outwards along the ray
(perhaps using the newton's method of mandel-exray.pdf in reverse),
repeat until no more bits left.  you move C a fractional iteration count
each time, and collect bits when crossing integer dwell boundaries
[Claude Heiland-Allen](http://mathr.co.uk/blog/)
>>>>











# See also
* [Computing external dynamic angle](https://gitlab.com/adammajewski/dynamic_external_angle)
* [parameter ray_in using Newton method and mpfr library in c ](https://gitlab.com/c_files/parameter_ray_in_newton_mpfr)
* [NonInteractive Parameter Ray_In MPFR](https://gitlab.com/adammajewski/NonInteractiveParameterRayInMPFR)


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

