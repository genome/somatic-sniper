#include "sniper_util.h"
#include <math.h>
#include <stdlib.h>

//phred space log add
double logAdd(double a, double c) {
    if(a < c) {
        return (a + logPhred(1.0 + expPhred(c - a) ) );
    } 
    else {
        return (c + logPhred(1.0 + expPhred(a - c) ) );
    }
}

double logSubtract(double a, double c) {
    if(a < c) {
        return (a + logPhred(1.0 - expPhred(c - a) ) );
    } 
    else {
        return (c + logPhred(1.0 - expPhred(a - c) ) );
    }

}
