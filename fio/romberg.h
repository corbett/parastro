/*
 ** Romberg integrator for an open interval.
 */
double dRombergO(void *CTX, double (*func)(void *, double), double a,
                 double b, double eps);
