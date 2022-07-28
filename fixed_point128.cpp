// fixed_point128.cpp : This file contains the 'main' function. Program execution begins and ends there.
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
    f3 = fp128::abs(-f3);
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
    fixed_point128<4> f1 = 1.0;
    fixed_point128<4> f2 = 1.4142135623730951; // == sqrt of 2
    fixed_point128<4> f1_sq = f1 * f1;
    fixed_point128<4> f2_sq = f2 * f2;
    printf("f1: %0.15lf\n", (double)f1);
    printf("f1 * f1: %0.15lf\n", (double)f1_sq);
    printf("f2: %0.15lf\n", (double)f2);
    printf("f2 * f2: %0.15lf\n", (double)f2_sq);
    f1 = 0.00001;
    printf("f1: %0.15lf\n", (double)f1);
    f1 *= 100;
    printf("f1 *= 100: %0.15lf\n", (double)f1);
    f1 = f1 * 100;
}

void test_division()
{
    printf("\nTest Division\n");
    fixed_point128<10> f1 = 1.0 / 3.0;
    printf("f1: %0.15lf\n", (double)f1);
    fixed_point128<10> f2 = f1 / 2.0;
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

}

void test_precision()
{
    printf("\nTest Precision\n");
    fixed_point128<10> f = 1;
    for (int i = 0; i < 118; ++i)         {
        f /= 2;
        printf("exp = -%d, f: %0.38lf\n", i, (double)f);
    }
}


void test_string()
{
    printf("\nTest string\n");
    fixed_point128<5> f1 = 1.0 / (3.0 * (1 << 21));
    std::string s = f1;
    printf("f1: %s\n", s.c_str());
}
void test_functions()
{
    printf("\nTest functions\n");
    double values[] = {3.45, -7.5};
    for (int i = 0; i < 2; ++i) {
        fixed_point128<5> f1 = values[i];
        printf("f1: %0.15lf\n", (double)f1);
        printf("floor(f1): %0.15lf\n", (double)floor(f1));
        printf("ciel(f1): %0.15lf\n", (double)ciel(f1));
        printf("abs(f1): %0.15lf\n", (double)abs(f1));
    }
}


int main()
{
    test_conversion();
    test_addition();
    test_shift();
    test_multiplication();
    test_division();
    test_precision();
    test_functions();
    test_string();
}

