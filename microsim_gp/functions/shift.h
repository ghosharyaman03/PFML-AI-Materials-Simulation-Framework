#ifndef APPLY_SHIFTY_H
#define APPLY_SHIFTY_H

void apply_shiftY(struct fields* gridinfo_w, long INTERFACE_POS_GLOBAL) {
  //Shift by one cell in the negative y-direction
  long x, y, z, dim;
  long gidy;
  double chemical_potential;
  double c[NUMCOMPONENTS-1];
  
  for(x=0; x < workers_mpi.rows_x; x++) {
    for(z=0; z < workers_mpi.rows_z; z++) {
      for (y=0; y <= (workers_mpi.rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y++) {
        gidy = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
        for (b=0; b < NUMPHASES; b++) {
          gridinfo_w[gidy].phia[b] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].phia[b];
        }
        if (FUNCTION_F != 5) {
          for (k=0; k < NUMCOMPONENTS-1; k++) {
            gridinfo_w[gidy].compi[k] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].compi[k];
          }
          for (k=0; k < NUMCOMPONENTS-1; k++) {
            gridinfo_w[gidy].composition[k] = gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].composition[k];
          }
        }
        gridinfo_w[gidy].temperature = gridinfo_w[gidy + (INTERFACE_POS_GLOBAL-shiftj)].temperature;
        if (LBM) {
          for (k=0;k<SIZE_LBM_FIELDS -4 ; k++) {
            lbm_gridinfo_w[gidy].fdis[k] = lbm_gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].fdis[k];
          }
          lbm_gridinfo_w[gidy].rho = lbm_gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].rho;
          for (dim=0; dim < DIMENSION; dim++) {
            lbm_gridinfo_w[gidy].u[dim] = lbm_gridinfo_w[gidy+(INTERFACE_POS_GLOBAL-shiftj)].u[dim];
          }
        }
      }

      if (FUNCTION_F != 5) {
        if (workers_mpi.lasty==1) {
          for (y=(workers_mpi.rows_y-4-(INTERFACE_POS_GLOBAL-shiftj)); y<=(workers_mpi.rows_y-4); y++) {

            gidy = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            for (b=0; b < NUMPHASES-1; b++) {
              gridinfo_w[gidy].phia[b] = 0.0;
            }
            gridinfo_w[gidy].phia[NUMPHASES-1] = 1.0;
            gridinfo_w[gidy].temperature = gridinfo_w[gidy-1].temperature + GRADIENT; 

            // if(!EIDT_SWITCH){
              for (k=0; k < NUMCOMPONENTS-1; k++) {
                c[k] = cfill[NUMPHASES-1][NUMPHASES-1][k];
              }
              Mu(c, Tfill, NUMPHASES-1, gridinfo_w[gidy].compi);
            // }
          }

          // if(EIDT_SWITCH){
          //   long y_start = workers_mpi.rows_y-4-(INTERFACE_POS_GLOBAL-shiftj) ;
          //   gidy         = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_start;

          //   for (y=y_start+1; y<=workers_mpi.rows_y-1; y++){
          //     long gidy_to = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

          //     for (k=0; k< NUMCOMPONENTS-1; k++) {
          //       EIDT.U0 = gridinfo_w[gidy].composition[k];    
          //       EIDT.B  = (EIDT.interface_velocity/Diffusivity[NUMPHASES-1][k][k]) * log( EIDT.U0 / EIDT.comp_far_field[k] ) ;
            
          //       gridinfo_w[gidy_to].composition[k] = EIDT.U0 * exp( -EIDT.B *(y - y_start)*deltay);

          //     }  
          //     Mu(gridinfo_w[gidy_to].composition, Tfill, NUMPHASES-1, gridinfo_w[gidy_to].compi);     
          //   }
          // }

          /**
           * The new fluid entering the system is stationary
           * sanjeev solver name alpha_beta
           */
          for (y=(workers_mpi.rows_y-1-(INTERFACE_POS_GLOBAL-shiftj)); y<=(workers_mpi.rows_y-1); y++) {
            gidy = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            if (LBM) {
              for (k=0;k<SIZE_LBM_FIELDS -4 ; k++) {
                lbm_gridinfo_w[gidy].fdis[k] = 0.0;
              }
              lbm_gridinfo_w[gidy].rho = rho_LBM[NUMPHASES-1];
              for (dim=0; dim < DIMENSION; dim++) {
                lbm_gridinfo_w[gidy].u[dim] = 0.0;
              }
            }
          }
        }
      }
    }
  }
  if(taskid ==  MASTER){
  printf("\n- [%d] shifting done!", taskid) ;
  }
//   init_propertymatrices(T);
}
long check_SHIFT(long x) {
  long center;
  long y, z;
  long INTERFACE_POS_MAX = 0;
  for (z=0; z < workers_mpi.rows_z; z++) {
    for (y=1; y <=(workers_mpi.rows_y-1); y++) {
  //     center =  gidy   + (x)*numy[levels];
      center = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y; 
  //     printf("center=%ld\n",center);
      if ((gridinfo_w[center-1].phia[NUMPHASES-1]-(1.0-gridinfo_w[center-1].phia[NUMPHASES-1]) < 0.0) 
        && (gridinfo_w[center].phia[NUMPHASES-1]-(1.0-gridinfo_w[center].phia[NUMPHASES-1]) > 0.0) ) {
        if (y > INTERFACE_POS_MAX) {
          INTERFACE_POS_MAX = y;
        }
      }
    }
  }
  if (INTERFACE_POS_MAX > 0) {
    return INTERFACE_POS_MAX - workers_mpi.offset_y + workers_mpi.offset[Y];
  } else {
    return 0;
  }
}
#endif
