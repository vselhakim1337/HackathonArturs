#include <iostream>
#include <math.h>

double hellosin(double a)
{
    return sin(a);
}

int main()
{
    double a = M_PI/2.0;
    double res = hellosin(a);
    return 0;
}