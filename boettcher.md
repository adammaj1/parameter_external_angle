
# The Böttcher function
The Böttcher function:
* maps the complement of the Mandelbrot set (or the filled Julia set)  conformally to the complement of the closed unit disk.
* is a solution of [Boettcher's functional equation](https://en.wikipedia.org/wiki/B%C3%B6ttcher%27s_equation), so in case of [complex quadratic map](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial)

$`\Phi_c: \hat{\C} \setminus K_c &\to \hat{\C} \setminus \overline{\mathbb{D}}\`$

$`B(f(z)) = (B(z))^2`$

where:

 - B is Boettchers function
 - f is [complex quadratic map](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial)



# Notation
* $`\Phi = B`$ is a Bottcher map ( function)
* $`arg`$ is the angle ( argument, phase) of complex number
* $`arg_M`$ is exernal angle on the parameter plane
* $`arg_c`$ is the external angle on the dynamic plane
* $`\Phi_c`$ is [the Boettcher map ](https://en.wikipedia.org/wiki/External_ray#Dynamical_plane_.3D_z-plane) on the dynamic plane
* $`\Phi_M`$ is [the Boettcher map ](https://en.wikipedia.org/wiki/External_ray#Dynamical_plane_.3D_z-plane) on the parameter plane


# solutions of the Boettchers equation

Solutions: 
* closed form 
* aproximation ( series )

## closed form solutions

## aproximations




### asymptotic series approximation for the Böttcher function


[On the parameter plane](http://reference.wolfram.com/language/ref/MandelbrotSetBoettcher.html) : 

$`B_M(c) = c + \frac{1}{2} - \frac{1}{8c} + \frac{1}{16c^2} - \frac{5}{128c^3} + O(\frac{1}{c^4})`$

[On the dynamic plane](http://reference.wolfram.com/language/ref/JuliaSetBoettcher.html) : 

$`B_c(z) = z + \frac{0.05}{z} + \frac{0.02375}{z^3} - \frac{0.0036875}{z^5} + O(\frac{1}{z^6})`$


Different values are described by <ref>[https://books.google.pl/books?id=SvT_AwAAQBAJ&pg=PA49&lpg=PA49&dq=%22boettcher+function%22&source=bl&ots=KIAZpgX-9y&sig=r1OztpQT7ITgGWSwtBBG1ipvVyY&hl=pl&sa=X&ved=0ahUKEwj49_f2teDWAhWGa1AKHWu1BAkQ6AEIQjAD#v=onepage&q=%22boettcher%20function%22&f=false ]</ref> [Lauwerier, H. A. (Hendrik Adolf) in the book: Chaos by Arun V. Holden Princeton University Press, 14 lip 2014, page 79 ]()

### computing argument and radius separtely



>" since the argument of the Boettcher function is not a main value but a value obtained by retaining multivalency, in fact, numerical calculation using this limit expression formula is very difficult. However, an infinite series easily obtained by transforming the limit expression formula"  


Souichiro-Ikebe ( automatic translation)  



$` \Phi_c(z) = \lim_{n\to \infty} (f_c^n(z))^{2^{-n}}  = R e^{i \theta}  = e^U e^{i \theta}`$


[argument of Boettcher coordinate ( = external angle)](http://math-functions-1.watson.jp/sub1_spec_390.html#section060): 

$`\theta = \theta_c(z) = arg(\Phi_c(z)) = arg_c(z) = arg(z) + \sum_{n=1}^\infty \left( \frac{1}{2^n}*arg \left(1 + \frac{c}{f_c^{n-1}(z)^2}     \right ) \right )  `$


Potential ( absolute value = radius = magnitude):

$`U = U_c(z) = log|\Phi_c(z)| =  log|z| + \sum_{n=1}^\infty \left( \frac{1}{2^n}*log \left|1 + \frac{c}{f_c^{n-1}(z)^2}     \right | \right )  `$


Links
* [Algorithm](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/boettcher#ArgPhi_-_External_angle_-_angular_component_of_complex_potential)
* Böttcher function by Souichiro-Ikebe
  * [description](http://math-functions-1.watson.jp/sub1_spec_390.html#section060)
  * [How to process branch cut lines of Böttcher function](http://math-functions-1.watson.jp/sub4_math_020.html#section030)
  * [FractalRelated.m ](the Mathematica-Package file " FractalRelated.m " in "Code of special function graph" of "File created in Version 8 ") - package of the Mathematica code


### n

(?????)

From the definitions:  
$`\Phi_M(c) = \Phi_c(c)`$  
$`arg_M(c)  = arg(\Phi_M(c)) `$  
$`arg_c(z) = arg(\Phi_c(z)) `$  
so [external angle of point c on the parameter plane is equal to external angle of the point z=c on the dynamic plane](https://en.wikipedia.org/wiki/External_ray)  
$`arg_M(c) = arg_c(z=c) = arg(\Phi_c(z= c)) = arg(\Phi_M(c)) `$
[//]: # "\Phi_M(c) \overset{def}{=} \Phi_c(c)`"




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
 
