#ifndef GCVSPL_H
#define GCVSPL_H

#include "f2c.h"

int gcvspl_( doublereal *x, doublereal *y, integer *ny, doublereal *wx, doublereal *wy, integer *m, integer *n, integer *k, integer *md, doublereal *val, doublereal *c__, integer *nc, doublereal *wk, integer *ier );
doublereal splder_( integer *ider, integer *m, integer *n, doublereal *t, doublereal *x, doublereal *c__, integer *l, doublereal *q );

#endif /* GCVSPL_H */


