# What is the external angle ? 

[External angle](boettcher.md#computing-argument-and-radius-separtely):
* is a real number -[description by Robert P Munafo](https://mrob.com/pub/muency/externalangle.html) 
* it is the angle of [external ray](https://en.wikipedia.org/wiki/External_ray) which 
  * goes through exterior point c of the Mandelbrot set
  * lands on the boundary point c of the Mandelbrot set
* it is the angle of point from the exterior of Mandelbrot set. It is not arg(c) but angle of [the Boettcher coordinate](boettcher.md) of point  = $` arg_M(c) `$ 
* angular part of complex potential





# Names of the external angle
* external argument
* phase



# How external angle is measured ?


## Direction
The angle is measured in [the counterclockwise or anticlockwise direction](https://en.wikipedia.org/wiki/Clockwise) 

## Units and domain
* when measured in [turns](https://en.wikipedia.org/wiki/Turn_(geometry)) it is a real number between 0.0 and 1.0
* when measured in radians it is a real number between 0.0 and $`2\Pi`$
* when measured in degrees it is a real number between 0.0 and 360

## base and expansion of numerical value
* [nine-ways-to-display-a-floating-point-number](http://www.exploringbinary.com/nine-ways-to-display-a-floating-point-number/)
* [What Every Programmer Should Know About Floating-Point Arithmetic ](http://floating-point-gui.de/ What Every Programmer Should Know About Floating-Point Arithmetic )
* [Stackoverflow : Why Are Floating Point Numbers Inaccurate?](http://stackoverflow.com/questions/21895756/why-are-floating-point-numbers-inaccurate )
* [HOW TO WORK WITH ONE-DI MENSIONAL QUADRATIC MAPS  G. Pastor  , M. Romera, G. Álvarez, and F. Montoya ](http://www.iec.csic.es/~gerardo/publica/Pastor03.pdf) 
* [tut](http://kipirvine.com/asm/workbook/floating_tut.htm)
* [home school math  : The fascinating irrational numbers](http://www.homeschoolmath.net/teaching/irrational_numbers.php )

### decimal number (base = 10 ) 
* real number
  * ratio = fraction ( Finite Continued fraction )  = rational number ( if number can not be represented as a ratio then it is irrational number ) 
    * in lowest terms ( irreducible form ) : $`\tfrac{1}{21}`$
    * reducible form 
      * in explicit normalized form ( only when denominator is odd ): $`\tfrac{3}{63} = \tfrac{3}{2^{6}-1}`$
  * irrational number( infinite continued fraction )
  *  floating point number $`0.\overline{047619}`$
    * finite expansion
    * endless expansion 
      *  continue infinitely without repeating (in which case the number is called irrational = non-repeating non-terminating decimal numbers)
      * Recurring or repeating : (strictly) periodic ( preperiod = 0 , preiod > 0 )  or mixed = eventually periodic ( preperiod > 0 , period > 0 )   

### binary number ( base = 2 

* rational number ( ratio) $`\tfrac{1}{10101}`$
* real number
  * floating point number ( scientific notation )
    * Raw binary (  raw IEEE format ) 
    *  fixed point number ( notation) 
      * with repeating sequences : $`0.\overline{000011}`$
      * with endless expansion $`0.000011000011000011000011...`$





# technical notes
GitLab uses:
* the Redcarpet Ruby library for [Markdown processing](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md)
* [KaTeX](https://khan.github.io/KaTeX/) to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html). [Here is used version](https://github.com/gitlabhq/gitlabhq/blob/a0715f079c143a362a7f6157db45020b8432003e/vendor/assets/javascripts/katex.js)


## git ( gitlab)

```
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/parameter_external_angle.git
git add .
git commit -m "Initial commit"
git push -u origin master
```

