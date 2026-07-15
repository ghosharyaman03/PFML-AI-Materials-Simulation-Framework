#ifndef ANISO_GEN_H
#define ANISO_GEN_H

void anisotropy_42_dAdq (double *qab, double* dadq, long a, long b)\
{
  const double qx2 = qab[X]*qab[X];
  const double qx4 = qx2*qx2;
  
  const double qy2 = qab[Y]*qab[Y];
  const double qy4 = qy2*qy2;
  
  const double qz2 = qab[Z]*qab[Z];
  const double qz4 = qz2*qz2;
  
  const double q2  = qx2 + qy2 + qz2;
  const double q23 = q2*q2*q2;
  
  if (fabs(q2) > 1.0e-15) {
    dadq[X] = 16.0*dab[a][b]*qab[X]*(( - qx4 - qz4 + qx2*q2)/q23);
    dadq[Y] = 16.0*dab[a][b]*qab[Y]*(( + qx4 + qz4         )/q23);
    dadq[Z] = 16.0*dab[a][b]*qab[Z]*(( - qx4 - qz4 + qz2*q2)/q23);
  } else {
    dadq[X] = 0.0;
    dadq[Y] = 0.0;
    dadq[Z] = 0.0;
  }
}

double anisotropy_42_function_ac(double *qab, long a, long b)
{
  const double qx2 = qab[X]*qab[X];
  const double qx4 = qx2*qx2;
  
  const double qy2 = qab[Y]*qab[Y];
  const double qy4 = qy2*qy2;
  
  const double qz2 = qab[Z]*qab[Z];
  const double qz4 = qz2*qz2;
  
  const double q2  = qx2 + qy2 + qz2;

  double ac;
  if (fabs(q2) > 1.0e-15) {
   ac = 1.0 - dab[a][b]*(3.0 - 4.0*(qx4 + qz4)/(q2*q2));
  } else {
   ac = 1.0;
  }
  return ac;
}

void anisotropy_gen_dAdq(double *qab, double* dadq, long a, long b)
{
  if (fab[a][b] == 2.0)
  {
    anisotropy_02_dAdq(qab, dadq, a, b);
  }
  else if (fab[a][b] == 4.0)
  {
    anisotropy_01_dAdq(qab, dadq, a, b);
  }
  else if (fab[a][b] == 4.2)
  {
    anisotropy_42_dAdq(qab, dadq, a, b);
  }
  else
  {
    dadq[X] = 0.0;
    dadq[Y] = 0.0;
    dadq[Z] = 0.0;
  }
}

double anisotropy_gen_function_ac(double *qab, long a, long b)
{
  if (fab[a][b] == 2.0)
  {
    return anisotropy_02_function_ac(qab, a, b); 
  }
  else if (fab[a][b] == 4.0)
  {
    return anisotropy_01_function_ac(qab, a, b);
  }
  else if (fab[a][b] == 4.2)
  {
    return anisotropy_42_function_ac(qab, a, b);
  }
  else
  {
    return 1.0;
  }
}

#endif