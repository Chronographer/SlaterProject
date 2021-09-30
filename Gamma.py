"""Gamma function borrowed from http://en.literateprograms.org/Gamma_function_with_the_Lanczos_approximation_(Python)"""

from cmath import *

g = 7
lanczos_coef = [ \
     0.99999999999980993,
   676.5203681218851,
 -1259.1392167224028,
   771.32342877765313,
  -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
     9.9843695780195716e-6,
     1.5056327351493116e-7]

def gamma(z):
    z = complex(z)
    if z.real < 0.5:
        return pi / (sin(pi*z)*gamma(1-z))
    else:
        z -= 1
        x = lanczos_coef[0] + \
            sum(lanczos_coef[i]/(z+i)
                for i in range(1, g+2))
        t = z + g + 0.5
        return sqrt(2*pi) * t**(z+0.5) * exp(-t) * x

## GOOD! -- anything by Lanczos in numerical modeling is AOK.  Its like using
##          an algorithm by Donald Knuth in computer science.
## NOTE -- make gamma its own module (separate file).
## NOTE -- Hopefully you intend to test the above definition with the below
##         dicitonary!! (and include the test in the separate module.)
#Test - print out in an array the values of gamma(z) relevant to the
if __name__=="__main__":
    print("Hello world from the gamma function!\n")
    print("Testing the gamma function ...\n")
    listgamma = [gamma(1),gamma(2),gamma(3),gamma(3.7),gamma(4),gamma(4.2)]
    nstar = {1:1.0, 2:2.0, 3:3.0, 4:3.7, 5:4.0, 6:4.2}
    gammafunction = [1.0,1.0,2.0,4.1706517838,6.0,7.75668953579]
    listout = []
    for i in range(len(listgamma)):
        print(i, listgamma[i], gammafunction[i], gammafunction[i]-listgamma[i])
