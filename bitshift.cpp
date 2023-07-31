#include <iostream>

int main()
{
    double a = 256.0;
    //int64_t *aint = reinterpret_cast<int64_t*>(&a);
    int64_t *aint = (int64_t*)(&a);
    (*aint) >>= 1;
    double *b = (double*)(aint);
    std::cout << *aint << '\n';
    std::cout << *b << '\n';

    return 0;
}