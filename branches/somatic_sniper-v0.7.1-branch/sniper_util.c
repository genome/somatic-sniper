#include "sniper_util.h"
#include <math.h>

//phred space log add
int logAdd(int a, int c) {
    if(a < c) {
        return (a + logPhred(1.0 + expPhred(c - a) ) );
    } 
    else {
        return (c + logPhred(1.0 + expPhred(a - c) ) );
    }
}
