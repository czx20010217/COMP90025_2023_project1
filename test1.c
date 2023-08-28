// Online C compiler to run C program online
#include <stdio.h>

typedef struct {
    int par1, par2;
} params;

double
cell_cost (long int seed, params par)
{
    const unsigned long a = 16807;
    const unsigned long m = 2147483647;

    /* For debugging only */
    // return (seed);

    /* Real code */
    seed = -seed;       // Make high bits non-zero
    int res   = par.par1;
    int scale = par.par2;

    int cost;
    
    for (cost = 0; seed >> res != 0; cost++) {
        seed = (a * seed) % m;
    }

    return (10 + (cost >> (8 * sizeof(unsigned long) - res - scale))) / 10.0;
}

int main() {
    // Write C code here
    printf("Start\n");
    params par;
    par.par1 = 10;
    par.par2 = 5;
    printf("start compute\n");
    double cost = cell_cost(3.5, par);
    printf("end: %f\n", cost);

    return 0;
}