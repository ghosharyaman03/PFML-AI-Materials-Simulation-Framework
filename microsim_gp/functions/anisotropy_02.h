#ifndef ANISOTROPY_02_H_
#define ANISOTROPY_02_H_

void anisotropy_02_dAdq (double *qab, double* dadq, long a, long b) {
  double qx2 = qab[X]*qab[X];
  double qy2 = qab[Y]*qab[Y];
  double qz2 = qab[Z]*qab[Z];
  double q2  = qx2 + qy2 + qz2;
  double q22  = q2*q2;
  
  if (fabs(q2) > 1.0e-15) {
   dadq[X] = 2.0*dab[a][b] * qab[X] * (   -qy2 -2.0*qz2) / q22;
   dadq[Y] = 2.0*dab[a][b] * qab[Y] * (    qx2 +    qz2) / q22;
   dadq[Z] = 2.0*dab[a][b] * qab[Z] * (2.0*qx2 +    qy2) / q22;
  } else {
    dadq[X] = 0.0;
    dadq[Y] = 0.0;
    dadq[Z] = 0.0;
  }

}

double anisotropy_02_function_ac(double *qab, long a, long b) {
  double qx2 = qab[X]*qab[X];
  double qy2 = qab[Y]*qab[Y];
  double qz2 = qab[Z]*qab[Z];
  double q2  = qx2 + qy2 + qz2;
  double ac  = 0.0;
  
  if (fabs(q2) > 1.0e-15) {
   ac = 1.0 + dab[a][b]*(qz2 - qx2) / q2;
  } else {
   ac = 1.0;
  }
  return ac;
}



#endif
