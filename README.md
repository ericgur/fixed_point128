# The fixed_point128 Class Template
 A 128 bit fixed-point class template for fast, high precision calculations.
 I developed this code during my COVID sick leave in a possibly futile effort to combat long-COVID.
 This code is used in my [Mandelbrot just-for-fun project](https://github.com/ericgur/Mandelbrot). With double precision floats I could zoom the image to 2^44^, with the fixed-point, 2^113^ is now possible.   
 
 Modifying existing floating point code to use this class is a fairly low effort work.

## fixed_point128 Properties
 - Most operations are very fast. 1-10x slower than double precision. ~10x faster than MPIR at similar precision.
 - Up to 38 fraction digits (decimal) are supported.
 - Has a superset of integer and floating point functions including all standard C/C++ operators.
 - The single template paramter **\<I\>** allows the user to specify 1-64 bits for the integer part, the rest are allocated to the fraction.
 - An object can be created from all int/float types as well as from strings representing a float.
 - Supports converions from one template instance to another (2 instances with different **\<I\>** parameter).
 
 ## Dependencies and Perquisites
 - Visual Studio 2019+ (MSFT Compiler)
 - C++17 or newer
 - Standard C++ library.
 - 64 bit builds only
 - Windows only. *Maybe* possible to adapt to other build environments.
 
 ***Note:*** Not tested with other compilers.

 ### The following standard math functions are supported:
 - fabs
 - floor, ceil
 - trunc, round
 - fmod
 - modf
 - sqrt with selectable precision/performance
 - sin, cos, asin, acos
 - tan, atan
 - exp, exp2, expm1
 - pow
 - log, log2, log10, logb, log1p
  ### Additional functions:
 - lzcnt128 - left zero count
 - fact_reciprocal - factorial reciprocal
 - reciprocal

 ## Visual Studio Debugger Support
 Add the file **fixed_point128.natvis** to your Visual Studio project.
 Now when debugging you'll see a pretty print of the objects value (limited to double precision).

 ## Code Examples

 ### Mandelbrot main loop
 Modified code snippet from my Mandelbrot C++/MFC project: https://github.com/ericgur/Mandelbrot
 Allows plotting the Mandelbrot (or Julia set) with a zoom of 2^113^.

    fixed_point128<8> x, y; // coordinate arguments.
    fixed_point128<8> radius = 2, radius_sq = radius * radius;
    fixed_point128<8> usq = 0, vsq = 0, u = 0, v = 0, tmp, uv, modulus = 0;
    
    // Find how many iterations are needed to have the (x,y) coordinates diverge (absolute value > 2)
    while (modulus < radius_sq && ++iter < MAX_ITERATION) {
        // real
        tmp = usq - vsq + xc;

        // imaginary
        //v = 2.0 * (u * v) + y;
        v = ((u * v) << 1) + yc;
        u = tmp;
        usq = u * u;
        vsq = v * v;
        // check uv vector amplitude is smaller than 2
        modulus = usq + vsq;
    }

## Acknologements
The function div_32bit (multi-precision integer division) is derived from the book *"Hacker's Delight"* 2nd Edition by Henry S. Warren Jr. 
It was converted to 32 bit operations and mdified a bit. The algorithm is an implementation of Knuth's "Algorithm D" from the book *"The Art of Computer Pogramming"*.

The functions log, log2, log10 are derived from Dan Moulding's code: https://github.com/dmoulding/log2fix

The function sqrt is based on the book *"Math toolkit for real time programming"* by Jack W. Crenshaw.
