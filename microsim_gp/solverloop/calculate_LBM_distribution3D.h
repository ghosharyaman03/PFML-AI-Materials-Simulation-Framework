#ifndef LBM_3D_FUNCTIONS_SOLVERLOOP
#define LBM_3D_FUNCTIONS_SOLVERLOOP

double w_3D[27] = {0.2962962962962963, 0.07407407407407407, 0.07407407407407407, 0.07407407407407407, 0.07407407407407407, 0.07407407407407407, 0.07407407407407407, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.018518518518518517, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629, 0.004629629629629629};

double vgrid_3D[3][27] = { 
  {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1},
  {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 1,-1, 1,-1,-1, 1, 1,-1},
  {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1}};


void calculate_distribution_function_3D(long x) {
  
  long dim, k, l, a, y,z, tmpVar = 0;
  double feq = 0;
  double magu = 0, magncap;
  double c_u = 0;
  double GD[3] = {0.0,0.0,0.0};
  double GB[3] = {0.0,0.0,0.0};
  double GS[3] = {0.0,0.0,0.0};
  double G1[3] = {0.0,0.0,0.0};
  double Gfdis = 0.0;
  double gradphi[3] = {0.0,0.0,0.0};
  double n[3]= {0.0,0.0,0.0};
  double rho = 0 ;

  for(z=0; z<=workers_mpi.end[Z]+2; z++){
    for(y=0; y<=workers_mpi.end[Y]+2; y++){

      center     = y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
      front      = y + z*workers_mpi.rows_y + (x+1)*workers_mpi.layer_size;
      back       = y + z*workers_mpi.rows_y + (x-1)*workers_mpi.layer_size;
      right      = y+1 + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
      left       = y-1 + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
      top        = y + (z+1)*workers_mpi.rows_y + x*workers_mpi.layer_size;
      bottom     = y + (z-1)*workers_mpi.rows_y + x*workers_mpi.layer_size;

      if(x > 0 && x <= workers_mpi.end[X] +2){
        gradphi[X] = 0.5*( gridinfo_w[front].phia[NUMPHASES-1] - gridinfo_w[back].phia[NUMPHASES-1] )/deltax ;
      }

      if(y > 0 && y <= workers_mpi.end[Y] +2){
        gradphi[Y] = 0.5*( gridinfo_w[right].phia[NUMPHASES-1] - gridinfo_w[left].phia[NUMPHASES-1] )/deltay ;
      }

      if(z > 0 && z <= workers_mpi.end[Z] +2){
        gradphi[Z] = 0.5*( gridinfo_w[top].phia[NUMPHASES-1] - gridinfo_w[bottom].phia[NUMPHASES-1] )/deltaz ;
      }
      
    
      magncap =  gradphi[X]*gradphi[X] + gradphi[Y]*gradphi[Y] + gradphi[Z]*gradphi[Z] + 1e-6 ;
      magncap = 1.0 / magncap ;
      n[X] = gradphi[X]*( magncap );
      n[Y] = gradphi[Y]*( magncap );
      n[Z] = gradphi[Z]*( magncap );

      // n cap is calculated
      
      // calculate the magnitude of the u field
      magu = 0.0;
      for (dim=0; dim<DIMENSION; dim++) {
        magu += lbm_gridinfo_w[center].u[dim]*lbm_gridinfo_w[center].u[dim];
      }
      
      // calculate the force terms
      // Calculates the drag force
      tmpVar = -0.5/(W_0*W_0)*nu*10*(1.0-gridinfo_w[center].phia[NUMPHASES-1])*(1.0-gridinfo_w[center].phia[NUMPHASES-1]) ;
      for(dim=0; dim < DIMENSION; dim++) {
        // GD[dim] = -0.5/(W_0*W_0)*nu*10*(1.0-gridinfo_w[center].phia[NUMPHASES-1])*(1.0-gridinfo_w[center].phia[NUMPHASES-1])*lbm_gridinfo_w[center].u[dim];
        GD[dim] = tmpVar*lbm_gridinfo_w[center].u[dim];
        GB[dim] = 0.0;
        GS[dim] = 0.0;
      }
      
      // Calculates the buoyancy force along Y
      tmpVar = -0.50*lbm_gridinfo_w[center].rho*g_y*(gridinfo_w[center].phia[NUMPHASES-1]) ;
      for (k=0; k < NUMCOMPONENTS-1; k++) {
        GB[X]  = 0.0 ;
        // GB[Y] += -0.50*lbm_gridinfo_w[center].rho*g_y*beta_c[k]*(gridinfo_w[center].composition[k]-ceq[NUMPHASES-1][NUMPHASES-1][k]) * (gridinfo_w[center].phia[NUMPHASES-1]);
        GB[Y] += tmpVar*beta_c[k]*(gridinfo_w[center].composition[k]-ceq[NUMPHASES-1][NUMPHASES-1][k]) ;
        GB[Z]  = 0.0 ;
      }

      // Calculates the shrinkage force
      for (dim=0; dim < DIMENSION; dim++) {
        for (a=0; a < NUMPHASES-1; a++) {
          GS[dim] += -((rho_LBM[a]-rho_LBM[NUMPHASES-1])/rho_LBM[NUMPHASES-1])*n[dim]*(gridinfo_w[center].deltaphi[a])/deltat;
        }
      }

      // for each discrete velocity dimension
      //double sigma_fk = 0 ;
      for (k=0; k < SIZE_LBM_FIELDS-4; k++) {
        // set these to zero
        c_u  = 0.0, Gfdis = 0.0;

        // calculate c_u
        for (dim=0; dim < DIMENSION; dim++) {
          c_u += lbm_gridinfo_w[center].u[dim]*vgrid_3D[dim][k];
        }

        // calculate total G force
        for (dim=0; dim<DIMENSION; dim++) {
          G1[dim]= w_3D[k]*(3.0*(vgrid_3D[dim][k]-lbm_gridinfo_w[center].u[dim]) + 9*c_u*vgrid_3D[dim][k]);
          Gfdis += G1[dim]*( GD[dim]+GB[dim]+GS[dim] );
        }
        
        // Calculate f_eq 
        feq = lbm_gridinfo_w[center].rho*w_3D[k]*(1.0 + 3.0*c_u + 4.5*c_u*c_u - 1.50*magu);

        // Update f(x+u.dx, t+dt) = (1-W)*f(x,t) + W*f_eq + G_f*dt
        // omega = 1.0/(3.0*nu_lbm + 0.5); defined in initialize variables
        lbm_gridinfo_w[center].fdis[k] = (1.0-omega)*lbm_gridinfo_w[center].fdis[k] + omega*feq + Gfdis*dt;

        //  G1=%le; GD=%le; GB=%le; GS=%le
        
        //if( ((feq < 0 )|| (lbm_gridinfo_w[center].fdis[k] < 0) ) && !isboundary){
        //   printf("\n-  [%d] WW (%ld,%ld,%ld)[%ld] feq=:%le; f=%le; Gf=%le;",taskid, x , y, z, k, feq, lbm_gridinfo_w[center].fdis[k], Gfdis) ;
        //}
      }
    }
  }
  
}

void calculate_streaming_row_forward_3D(long x) {
  long y,z, center, back;
  
  for(z=1; z<=workers_mpi.end[Z]+2; z++){
  for(y=1; y<=workers_mpi.end[Y]+2; y++){
    center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
    back       =  y + z*workers_mpi.rows_y + (x-1)*workers_mpi.layer_size;

    lbm_gridinfo_w[back].fdis[2]   =  lbm_gridinfo_w[center].fdis[2];
    lbm_gridinfo_w[back].fdis[8]   =  lbm_gridinfo_w[center].fdis[8];
    lbm_gridinfo_w[back].fdis[10]  =  lbm_gridinfo_w[center].fdis[10];
    lbm_gridinfo_w[back].fdis[14]  =  lbm_gridinfo_w[center].fdis[14];
    lbm_gridinfo_w[back].fdis[16]  =  lbm_gridinfo_w[center].fdis[16];
    lbm_gridinfo_w[back].fdis[20]  =  lbm_gridinfo_w[center].fdis[20];
    lbm_gridinfo_w[back].fdis[22]  =  lbm_gridinfo_w[center].fdis[22];
    lbm_gridinfo_w[back].fdis[24]  =  lbm_gridinfo_w[center].fdis[24];
    lbm_gridinfo_w[back].fdis[25]  =  lbm_gridinfo_w[center].fdis[25];
  
  }
  }
  //printf("\n") ;
}

void calculate_streaming_row_back_3D(long x) {
  long y,z, center, back;
  
  for(z=1 ;z<=workers_mpi.end[Z]+2; z++){
  for(y=1; y<=workers_mpi.end[Y]+2; y++){
      center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
      back       =  y + z*workers_mpi.rows_y + (x-1)*workers_mpi.layer_size;

      //printf("\nRowBack(%ld,%ld,%ld)[%ld]->[%ld]",x, gidy/workers_mpi.rows_y, gidy%workers_mpi.rows_y, back, center);
      lbm_gridinfo_w[center].fdis[1]  =  lbm_gridinfo_w[back].fdis[1];
      lbm_gridinfo_w[center].fdis[7]  =  lbm_gridinfo_w[back].fdis[7];
      lbm_gridinfo_w[center].fdis[9]  =  lbm_gridinfo_w[back].fdis[9];
      lbm_gridinfo_w[center].fdis[13] =  lbm_gridinfo_w[back].fdis[13];
      lbm_gridinfo_w[center].fdis[15] =  lbm_gridinfo_w[back].fdis[15];
      lbm_gridinfo_w[center].fdis[19] =  lbm_gridinfo_w[back].fdis[19];
      lbm_gridinfo_w[center].fdis[21] =  lbm_gridinfo_w[back].fdis[21];
      lbm_gridinfo_w[center].fdis[23] =  lbm_gridinfo_w[back].fdis[23];
      lbm_gridinfo_w[center].fdis[26] =  lbm_gridinfo_w[back].fdis[26];
    }
  }
  //printf("\n") ;
}

void calculate_streaming_col_forward_3D(long x) {
  long y,z, center, left;
  
  for(z=1; z<=workers_mpi.end[Z]+2; z++){
  for(y=1; y<=workers_mpi.end[Y]+2; y++){
    center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
    left       =  y-1 + z*workers_mpi.rows_y + x*workers_mpi.layer_size;

    //printf("\nColFwd(%ld,%ld,%ld)[%ld]->[%ld]",x, gidy/workers_mpi.rows_y, gidy%workers_mpi.rows_y, center, left);
    lbm_gridinfo_w[left].fdis[4]   =  lbm_gridinfo_w[center].fdis[4];
    lbm_gridinfo_w[left].fdis[8]   =  lbm_gridinfo_w[center].fdis[8];
    lbm_gridinfo_w[left].fdis[12]  =  lbm_gridinfo_w[center].fdis[12];
    lbm_gridinfo_w[left].fdis[13]  =  lbm_gridinfo_w[center].fdis[13];
    lbm_gridinfo_w[left].fdis[18]  =  lbm_gridinfo_w[center].fdis[18];
    lbm_gridinfo_w[left].fdis[20]  =  lbm_gridinfo_w[center].fdis[20];
    lbm_gridinfo_w[left].fdis[22]  =  lbm_gridinfo_w[center].fdis[22];
    lbm_gridinfo_w[left].fdis[23]  =  lbm_gridinfo_w[center].fdis[23];
    lbm_gridinfo_w[left].fdis[26]  =  lbm_gridinfo_w[center].fdis[26];
  }}
  //printf("\n") ;
}

void calculate_streaming_col_back_3D(long x) {
  long y,z, center, left;
    
  for(z=1; z <=workers_mpi.end[Z]+2; z++){
  for(y=workers_mpi.end[Y]+2;   y > 1; y--){
    center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
    left       =  y-1 + z*workers_mpi.rows_y + x*workers_mpi.layer_size;


    //printf("\nColBack(%ld,%ld,%ld)[%ld]->[%ld]",x, gidy/workers_mpi.rows_y, gidy%workers_mpi.rows_y, left, center);
    lbm_gridinfo_w[center].fdis[3]  =  lbm_gridinfo_w[left].fdis[3];
    lbm_gridinfo_w[center].fdis[7]  =  lbm_gridinfo_w[left].fdis[7];
    lbm_gridinfo_w[center].fdis[11] =  lbm_gridinfo_w[left].fdis[11];
    lbm_gridinfo_w[center].fdis[14] =  lbm_gridinfo_w[left].fdis[14];
    lbm_gridinfo_w[center].fdis[17] =  lbm_gridinfo_w[left].fdis[17];
    lbm_gridinfo_w[center].fdis[19] =  lbm_gridinfo_w[left].fdis[19];
    lbm_gridinfo_w[center].fdis[21] =  lbm_gridinfo_w[left].fdis[21];
    lbm_gridinfo_w[center].fdis[24] =  lbm_gridinfo_w[left].fdis[24];
    lbm_gridinfo_w[center].fdis[25] =  lbm_gridinfo_w[left].fdis[25];
  }}
  //printf("\n") ;
}

void calculate_streaming_Z_forward_3D(long x) {
  long y,z, center, bottom;
  
  for(z=1; z<=workers_mpi.end[Z]+2; z++){
  for(y=1; y<=workers_mpi.end[Y]+2; y++){
    center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
    bottom     =  y + (z-1)*workers_mpi.rows_y + x*workers_mpi.layer_size;

    lbm_gridinfo_w[bottom].fdis[6]  =  lbm_gridinfo_w[center].fdis[6];
    lbm_gridinfo_w[bottom].fdis[10] =  lbm_gridinfo_w[center].fdis[10];
    lbm_gridinfo_w[bottom].fdis[12] =  lbm_gridinfo_w[center].fdis[12];
    lbm_gridinfo_w[bottom].fdis[15] =  lbm_gridinfo_w[center].fdis[15];
    lbm_gridinfo_w[bottom].fdis[17] =  lbm_gridinfo_w[center].fdis[17];
    lbm_gridinfo_w[bottom].fdis[20] =  lbm_gridinfo_w[center].fdis[20];
    lbm_gridinfo_w[bottom].fdis[21] =  lbm_gridinfo_w[center].fdis[21];
    lbm_gridinfo_w[bottom].fdis[24] =  lbm_gridinfo_w[center].fdis[24];
    lbm_gridinfo_w[bottom].fdis[26] =  lbm_gridinfo_w[center].fdis[26];
  }}
  //printf("\n") ;
}


void calculate_streaming_Z_back_3D(long x) {
  long y,z, center, bottom;
  
  for(z=workers_mpi.end[Z]+2; z>1; z--){
  for(y=1; y<=workers_mpi.end[Y]+2; y++){
    center     =  y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;
    bottom     =  y + (z-1)*workers_mpi.rows_y + x*workers_mpi.layer_size;

    //printf("\nZBack(%ld,%ld,%ld)[%ld]->[%ld]",x, gidy/workers_mpi.rows_y, gidy%workers_mpi.rows_y, bottom, center);
    lbm_gridinfo_w[center].fdis[5]  =  lbm_gridinfo_w[bottom].fdis[5];
    lbm_gridinfo_w[center].fdis[9]  =  lbm_gridinfo_w[bottom].fdis[9];
    lbm_gridinfo_w[center].fdis[11] =  lbm_gridinfo_w[bottom].fdis[11];
    lbm_gridinfo_w[center].fdis[16] =  lbm_gridinfo_w[bottom].fdis[16];
    lbm_gridinfo_w[center].fdis[18] =  lbm_gridinfo_w[bottom].fdis[18];
    lbm_gridinfo_w[center].fdis[19] =  lbm_gridinfo_w[bottom].fdis[19];
    lbm_gridinfo_w[center].fdis[22] =  lbm_gridinfo_w[bottom].fdis[22];
    lbm_gridinfo_w[center].fdis[23] =  lbm_gridinfo_w[bottom].fdis[23];
    lbm_gridinfo_w[center].fdis[25] =  lbm_gridinfo_w[bottom].fdis[25];
  }}
  //printf("\n") ;
}



void compute_field_variables_LBM_3D(long x) {
  
  long k,dim,gidy;
  double rho, u[3];
  
  // for(z=0; z<=workers_mpi.end[Z]+3; z++){
  // for(y=0; y<=workers_mpi.end[Y]+3; y++){

  // center     = y + z*workers_mpi.rows_y + x*workers_mpi.layer_size;

  for (gidy=1; gidy < workers_mpi.layer_size; gidy++){
    
    center        =  gidy   + (x)*workers_mpi.layer_size;

    rho = 0.0;
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;

    for (k=0; k < SIZE_LBM_FIELDS-4 ; k++) {
      rho   += lbm_gridinfo_w[center].fdis[k];
      for(dim=0; dim<DIMENSION;dim++){
        u[dim]  += lbm_gridinfo_w[center].fdis[k]*vgrid_3D[dim][k];
      }
    }

    // calculate rho and u field
    lbm_gridinfo_w[center].rho = rho;

    
    if (rho > workers_max_min.rho_max) {
      workers_max_min.rho_max = rho;
    } else if(rho < workers_max_min.rho_min) {
      workers_max_min.rho_min = rho;
    }

    for(dim=0; dim<DIMENSION;dim++){
      lbm_gridinfo_w[center].u[dim] = u[dim]/rho;
    }

  }

  //}

}



#endif
