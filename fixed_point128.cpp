/***********************************************************************************
    MIT License

    Copyright (c) 2022 Eric Gur (ericgur@iname.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
************************************************************************************/

// fixed_point128.cpp : test executable for the fixed_point128 class template
//
#include <stdio.h>
#include "fixed_point128.h" 

using namespace fp128;

void test_conversion()
{
    printf("\nTest Conversion\n");
    double a1 = 5.5;
    double a2 = 1.0 / (1ull << 22);
    double a3 = 3.000003;
    fixed_point128<4> f1 = a1;
    fixed_point128<4> f2 = a2;
    fixed_point128<32> f3 = a3;
    double r1 = f1;
    double r2 = f2;
    double r3 = f3;
    
    printf("a1: %0.15lf --> %0.15lf\n", a1, r1);
    printf("a2: %0.15lf --> %0.15lf\n", a2, r2);
    printf("a3: %0.15lf --> %0.15lf\n", a3, r3);
    f3 = fp128::fabs(-f3);
}

void test_addition()
{
    printf("\nTest Addition\n");
    fixed_point128<4> f1 = 5.123456789012345;
    fixed_point128<4> f2 =  7.0 - 5.123456789012345;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    ++f1;
    printf("++f1 --> %0.15lf\n", (double)f1);
    f1++;
    printf("f1++ --> %0.15lf\n", (double)f1);
    fixed_point128<4> f3 = f1 + f2;
    printf("f1+f2 --> %0.15lf\n", (double)f3);
    f3 = f2 + f1;
    printf("f2+f1 --> %0.15lf\n", (double)f3);
    f1 = -f1;
    printf("f1 = -f1 --> %0.15lf\n", (double)f1);
    f3 = f1 + f2;
    printf("f1+f2 --> %0.15lf\n", (double)f3);
    f3 = f1 - f2;
    printf("f1-f2 --> %0.15lf\n", (double)f3);
    f1 = 2.25; f2 = 4.54;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    f3 = f1 - f2;
    printf("f1-f2 --> %0.15lf\n", (double)f3);
    printf("Pi (double): %0.40lf\n", 3.14159265358979323846264338327950288419716939937510);
    printf("f1 = \"3.14159265358979323846264338327950288419716939937510\"\n"),
    f1 = "3.14159265358979323846264338327950288419716939937510"; // 50 first digits of pi
    std::string s = f1;
    printf("f1: %s\n", s.c_str());
}

void test_shift()
{
    printf("\nTest Shift\n");
    fixed_point128<4> f1 = 1.0;
    fixed_point128<4> f2 = 5.123456789012345;
    printf("f1: %0.15lf\n", (double)f1);
    f1 >>= 10;
    printf("f1 >>= 10 --> %0.15lf\n", (double)f1);
    f1 <<= 11;
    printf("f1 <<= 11 --> %0.15lf\n", (double)f1);
    printf("f2: %0.15lf\n", (double)f2);
    f2 >>= 65;
    printf("f2 >>= 65 --> %0.15lf\n", (double)f2);
    f2 <<= 66;
    printf("f2 <<= 66 --> %0.15lf\n", (double)f2);
}

void test_multiplication()
{
    printf("\nTest Multiplication\n");
    fixed_point128<16> f1 = 1.0;
    fixed_point128<16> f2 = 1.4142135623730951; // == sqrt of 2
    fixed_point128<16> f1_sq = f1 * f1;
    fixed_point128<16> f2_sq = f2 * f2;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f1 * f1: %0.15lf\n", (double)f1_sq);
    printf("f2: %0.15lf\n", (double)f2);
    printf("f2 * f2: %0.15lf\n", (double)f2_sq);
    f1 = 0.00001;
    printf("f1: %0.15lf\n", (double)f1);
    f1 *= -100;
    printf("f1 *= -100: %0.15lf\n", (double)f1);
    f1 = f1 * 100;
}

void test_division()
{
    printf("\nTest Division\n");
    fixed_point128<20> f1 = 1.0 / 3.0;
    printf("f1: %0.15lf\n", (double)f1);
    fixed_point128<20> f2 = f1 / 2.0;
    printf("f2 = f1 / 2.0: %0.15lf\n", (double)f2);
    f2 /= 32.0;
    printf("f2 = f1 / 32.0: %0.15lf\n", (double)f2);
    f1 = 0.01 / 3.0;
    printf("f1 = 0.01 / 3.0 %0.15lf\n", (double)f1);
    f2 = 1.0;
    printf("f2: %0.15lf\n", (double)f2);
    f2 /= f1;
    printf("f2 /= f1: %0.15lf\n", (double)f2);

    f1 = 8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 5;
    printf("f1 %%= 5: %0.15lf\n", (double)f1);
    f1 = 8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= -5;
    printf("f1 %%= -5: %0.15lf\n", (double)f1);
    f1 = -8;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 5;
    printf("f1 %%= 5: %0.15lf\n", (double)f1);
    f1 = 11;
    printf("f1: %0.15lf\n", (double)f1);
    f1 %= 0.7;
    printf("f1 %%= 0.7: %0.15lf\n", (double)f1);
    f1 = 11; 
    f1 %= -0.7;
    printf("f1 %%= -0.7: %0.15lf\n", (double)f1);
}

void test_precision()
{
    printf("\nTest Precision\n");
    fixed_point128<10> f = 1;
    for (int i = 1; i <= 118; ++i)         {
        f /= 2;
        printf("exp = -%d, f: %0.40lf\n", i, (double)f);
    }
}


void test_string()
{
    printf("\nTest string\n");
    double val = 1.0 / (3.0 * (1 << 21));
    fixed_point128<5> f1 = val;
    std::string s = f1;
    printf("double val:\n%0.80lf\n", val);
    printf("f1:\n%s\n", s.c_str());
    val = -5.12345678901234567890;
    f1 = val;
    s = f1;
    printf("double val:\n%0.80lf\n", val);
    printf("f1:\n%s\n", s.c_str());
}


void test_functions()
{
    printf("\nTest functions\n");
    double values[] = {3.45, -7.5};
    for (int i = 0; i < 2; ++i) {
        fixed_point128<10> f1 = values[i];
        printf("f1: %0.15lf\n", (double)f1);
        printf("floor(f1): %0.15lf\n", (double)floor(f1));
        printf("ciel(f1): %0.15lf\n", (double)ciel(f1));
        printf("fabs(f1): %0.15lf\n", (double)fabs(f1));
        fixed_point128<10> int_part;
        fixed_point128<10> f2 = modf(f1, &int_part);
        printf("modf(f1): int: %lf fraction %0.15lf\n", (double)int_part, (double)f2);
       
        f2 = 0.25;
        printf("f2: %0.15lf\n", (double)f2);
        printf("fmod(f1, f2): %0.15lf\n", (double)fmod(f1, f2));
        printf("\n");
    }
}


int main()
{
    test_precision();
    test_conversion();
    test_addition();
    test_shift();
    test_multiplication();
    test_division();
    test_functions();
    test_string();
}

