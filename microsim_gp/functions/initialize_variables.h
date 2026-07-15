#ifndef INITIALIAZE_VARIABLES_H_
#define INITIALIAZE_VARIABLES_H_

void initialize_variables()
{
  long a, i, j;
  for (a = 0; a < NUMPHASES; a++) {
    for (i = 0; i < NUMCOMPONENTS-1; i++) {
      for (j = 0; j < NUMCOMPONENTS-1; j++) {
        if (i==j) {
          muc[a][i][j]=2.0*A[a][i][j];
        } else {
          muc[a][i][j]=A[a][i][j];
        }
      }
    }
    matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
  }
  
  for (a=0;a<NUMPHASES;a++) {
    for (i=0; i < NUMCOMPONENTS-1; i++) {      
      dcbdT_phase[a][i] = 0.0;
      for (j=0; j < NUMCOMPONENTS-1; j++) {
        dcbdT_phase[a][i] += cmu[a][i][j]*(-dBbdT[a][j]);
      }
    }
  }

  if (USE_NEW_MALLOC)
  {    
    dcdmu       = ms_Malloc_2D((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
    dcdmu_phase = ms_Malloc_3D(NUMPHASES, NUMCOMPONENTS-1,NUMCOMPONENTS-1);
    Ddcdmu      = ms_Malloc_2D(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
    inv_dcdmu   = ms_Malloc_2D(NUMCOMPONENTS-1,NUMCOMPONENTS-1);

    deltamu     = ms_Malloc_1D(NUMCOMPONENTS-1);
    deltac      = ms_Malloc_1D(NUMCOMPONENTS-1);
    sum         = ms_Malloc_1D(NUMCOMPONENTS-1);
    divphi      = ms_Malloc_1D(NUMPHASES);
    lambda_phi  = ms_Malloc_1D(NUMPHASES);
    divflux     = ms_Malloc_1D(NUMCOMPONENTS-1);
    c_old       = ms_Malloc_1D(NUMCOMPONENTS-1);
    c_new       = ms_Malloc_1D(NUMCOMPONENTS-1);
    c_tdt       = ms_Malloc_1D(NUMCOMPONENTS-1);
    divjat      = ms_Malloc_1D(NUMCOMPONENTS-1);
  }
  else
  {
    dcdmu       = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
    dcdmu_phase = Malloc3M(NUMPHASES, NUMCOMPONENTS-1,NUMCOMPONENTS-1);
    Ddcdmu      = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
    inv_dcdmu   = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
    deltamu     = MallocV(NUMCOMPONENTS-1);
    deltac      = MallocV(NUMCOMPONENTS-1);
    sum         = MallocV(NUMCOMPONENTS-1);
    divphi      = MallocV(NUMPHASES);
    lambda_phi  = MallocV(NUMPHASES);
    divflux     = MallocV(NUMCOMPONENTS-1);
    c_old       = MallocV(NUMCOMPONENTS-1);
    c_new       = MallocV(NUMCOMPONENTS-1);
    c_tdt       = MallocV(NUMCOMPONENTS-1);
    divjat      = MallocV(NUMCOMPONENTS-1);
  }
  
  gridinfo_instance = (struct fields *)malloc(sizeof(*gridinfo_instance));
  allocate_memory_fields(gridinfo_instance);
  
  if (LBM) {
    if (DIMENSION ==2) {
      SIZE_LBM_FIELDS = 13;
    }
    if (DIMENSION == 3) {
      SIZE_LBM_FIELDS = 31;
    }
    lbm_gridinfo_w_instance = (struct lbm_fields *)malloc(sizeof(struct lbm_fields));
    allocate_memory_lbm_fields(lbm_gridinfo_w_instance);
    
    omega = 1.0/(3.0*nu_lbm + 0.5);

    if(taskid == MASTER){
      printf("\n- [0] LBM nu      = %12.12f", nu_lbm) ;
      printf("\n- [0] LBM omega   = %12.12f", omega) ;
      printf("\n- [0] LBM Tau     = %12.12f", 1.0/omega) ;
      printf("\n- [0] LBM Tau/dt  = %12.12f", (1.0/omega)/dt) ;
      printf("\n- [0] LBM |u|_max = %12.12f", 0.577*deltax/dt) ;
    }

  }

  EIDT.shift_interval       = 1;
  EIDT.interface_velocity   = 0;
  EIDT.B  = 0;
  EIDT.U0 = 1;
}

#endif
