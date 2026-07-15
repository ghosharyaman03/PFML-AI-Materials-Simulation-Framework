#ifndef BOUNDARY_MPI_H_
#define BOUNDARY_MPI_H_

void copyXZ(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo_w, char *field_type);
void copyXY(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo_w, char *field_type);
void copyYZ(struct bc_scalars *boundary, struct fields* gridinfo_w, char *field_type);
void constraint_uedge(struct bc_scalars *boundary, long y_this, long x_start, long x_end, long copy_to, long copy_from, long sign);

void constraint_uedge(struct bc_scalars *boundary, long y_edge, long x_start, long x_end, long copy_to, long copy_from, long sign) {
  long gidy_from;
  long gidy_to;
  long x,y,z;
  long j,m;
  long front, back, top, bottom, right, left;
  long copy_to_1, copy_from_1;
  struct symmetric_tensor eigen_strain;
  struct symmetric_tensor strain_grid;
  for (x=0; x < workers_mpi.rows_x; x++) {
    for (z=0; z < workers_mpi.rows_z; z++) {
      gidy_from         = x*workers_mpi.layer_size + z*workers_mpi.rows_y         + copy_from;
      gidy_to           = x*workers_mpi.layer_size + z*workers_mpi.rows_y         + copy_to;
      
      gridy_elast_first = x*workers_mpi.layer_size + z*workers_mpi.rows_y         + y_edge;
      
      if (x < workers_mpi.rows_x-1) {
        front           = (x+1)*workers_mpi.layer_size + z*workers_mpi.rows_y     + y_edge;
      } else {
        front           = (x)*workers_mpi.layer_size + z*workers_mpi.rows_y       + y_edge;
      }
      
      if (x > 0) {
        back            = (x-1)*workers_mpi.layer_size + z*workers_mpi.rows_y     + y_edge;
      } else {
        back            = x*workers_mpi.layer_size + z*workers_mpi.rows_y         + y_edge;
      }
      
      right             =  x*workers_mpi.layer_size + z*workers_mpi.rows_y        + y_edge + 1;
      left              =  x*workers_mpi.layer_size + z*workers_mpi.rows_y        + y_edge - 1;
      
      if (z < (workers_mpi.rows_z-1)) {
        top             =  x*workers_mpi.layer_size + (z+1)*workers_mpi.rows_y    + y_edge;
      } else {
        top             =  x*workers_mpi.layer_size + z*workers_mpi.rows_y        + y_edge;
      }
      
      if (z > 0) {
        bottom          =  x*workers_mpi.layer_size + (z-1)*workers_mpi.rows_y    + y_edge;
      } else {
        bottom          =  x*workers_mpi.layer_size + z*workers_mpi.rows_y        + y_edge;
      }
//       gridy_elast_next  = (x+1)*workers_mpi.layer_size + z*workers_mpi.rows_y + y_edge;
      stiffness         = calculate_stiffness_n(gridy_elast_first);
      eigen_strain      = calculate_eigen_strain(gridy_elast_first);
      
      
//       strain            = boundary[3].value[0]*stiffness.C11/(stiffness.C11*stiffness.C11 - stiffness.C12*stiffness.C12);
      strain_grid.xx    = ((0.5)*(iter_gridinfo_w[front].disp[X][2] - iter_gridinfo_w[back].disp[X][2]))
                                 - eigen_strain.xx;
//       strain_grid.xy    = (0.5)*(iter_gridinfo_w[front].disp[Y][2] - iter_gridinfo_w[back].disp[Y][2]);
    
                                 
      if (DIMENSION == 3) {
         strain_grid.zz   = (0.5)*((iter_gridinfo_w[top].disp[Z][2] - iter_gridinfo_w[bottom].disp[Z][2]))  
                                 - eigen_strain.zz;
                                 
//          strain_grid.yz   = 0.5*(iter_gridinfo_w[top].disp[Y][2] - iter_gridinfo_w[bottom].disp[Y][2]);
                                 
      } else {
        strain_grid.zz    = 0.0;
        strain_grid.yz    = 0.0; 
      }
                                 
      strain            = (boundary[3].value[0]/stiffness_phase[NUMPHASES-1].C44 - stiffness.C12*strain_grid.xx - stiffness.C12*strain_grid.zz + stiffness.C11*eigen_strain.yy)/(stiffness.C11);
      
      uy                = iter_gridinfo_w[gidy_from].disp[Y][2] + sign*2.0*strain;

      iter_gridinfo_w[gidy_to].disp[Y][2]           = uy;
      iter_gridinfo_w[gidy_to].disp[X][2]           = iter_gridinfo_w[gidy_from].disp[X][2];
      iter_gridinfo_w[gidy_to].disp[Z][2]           = iter_gridinfo_w[gidy_from].disp[Z][2];
      
      for(j=1; j<3; j++) {
        
        copy_to_1           = boundary[3].points[j];
        copy_from_1         = boundary[3].points[0];
        
        gidy_from           = x*workers_mpi.layer_size  +  z*workers_mpi.rows_y     + copy_from_1;
        gidy_to             = x*workers_mpi.layer_size  +  z*workers_mpi.rows_y     + copy_to_1;
        
        for (n=0; n < 3; n++) {
          iter_gridinfo_w[gidy_to].disp[n][2] = iter_gridinfo_w[gidy_from].disp[n][2];
        }
      }
    }
  }
}

void copyXZ(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo_w, char *field_type) {
  long gidy_from, gidy_to, y, a, k, y_start;
  long copy_from, copy_to;
  long j;
  long x, z;
  int m;
  if (strcmp(field_type, "PHI") == 0) {
   if(((boundary[0].type ==1) && (workers_mpi.firsty || workers_mpi.lasty)) || ((boundary[0].type == 3) && (workers_mpi.firsty && workers_mpi.lasty))) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[0].proxy[j];
        copy_to   = boundary[0].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to            = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
            for (a=0; a < NUMPHASES; a++) {
              gridinfo_w[gidy_to].phia[a]     = gridinfo_w[gidy_from].phia[a];
              gridinfo_w[gidy_to].deltaphi[a] = gridinfo_w[gidy_from].deltaphi[a];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "MU") == 0) {
   if(((boundary[1].type == 1) && (workers_mpi.firsty || workers_mpi.lasty)) || ((boundary[1].type == 3) && (workers_mpi.firsty && workers_mpi.lasty))) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[1].proxy[j];
        copy_to   = boundary[1].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to            = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
            for (k=0; k< NUMCOMPONENTS-1; k++) {
              gridinfo_w[gidy_to].compi[k]       = gridinfo_w[gidy_from].compi[k];
              gridinfo_w[gidy_to].composition[k] = gridinfo_w[gidy_from].composition[k];
            }
          }
        }
      }
    }
    if((boundary[1].type == 2) && (workers_mpi.firsty || workers_mpi.lasty)){
      if (workers_mpi.firsty) {
        y_first            = boundary[1].proxy[0] - 1;
        y_grid             = y_first;
      }
      if (workers_mpi.lasty) {
        y_last             = boundary[1].proxy[0] + 1;
        y_grid             = y_last;
      }
      for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_to   = boundary[1].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_to            = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
            gridy              = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_grid;

            Mu(ceq[NUMPHASES-1][NUMPHASES-1], Teq, NUMPHASES-1, gridinfo_w[gridy].compi);
            for (k=0; k< NUMCOMPONENTS-1; k++) {
              gridinfo_w[gidy_to].compi[k]       = gridinfo_w[gridy].compi[k];
              gridinfo_w[gidy_to].composition[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
            }            
          }
        }
      }
    }

    if(EIDT_SWITCH){
      if((boundary[1].type == 1) && (workers_mpi.lasty)){
        y_start = workers_mpi.rows_y - 6 ;
        for (y = y_start; y<=workers_mpi.rows_y-1; y++){
          for (x=x_start; x<=x_end; x++) {
            for (z=0; z < workers_mpi.rows_z; z++) {
              gidy    = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_start;  
              gidy_to = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

              for (k=0; k< NUMCOMPONENTS-1; k++) {
                
                EIDT.U0 = gridinfo_w[gidy].composition[k];    
                EIDT.B  = EIDT.interface_velocity / Diffusivity[NUMPHASES-1][k][k];

                gridinfo_w[gidy_to].composition[k] = EIDT.comp_far_field[k] - (EIDT.comp_far_field[k] - EIDT.U0) * exp( -EIDT.B *(y - y_start)*deltay);
              }

              Mu(gridinfo_w[gidy_to].composition, Tfill, NUMPHASES-1, gridinfo_w[gidy_to].compi);     
            }
          }
        }
      }
    }
  }

  if (strcmp(field_type, "T") == 0) {
   if(((boundary[2].type == 1) && (workers_mpi.firsty || workers_mpi.lasty)) || ((boundary[2].type == 3) && (workers_mpi.firsty && workers_mpi.lasty))) {
     for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[2].proxy[j];
        copy_to   = boundary[2].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from                       = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to                         = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
            gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
          }
        }
      }
    }
  }
  if (strcmp(field_type, "U") == 0) {
   if(((boundary[3].type == 1) && (workers_mpi.firsty || workers_mpi.lasty)) || ((boundary[3].type == 3) && (workers_mpi.firsty && workers_mpi.lasty)) ||((boundary[3].type == 2) && (workers_mpi.firsty && workers_mpi.lasty))) {
     for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[3].proxy[j];
        copy_to   = boundary[3].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from                       = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to                         = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
//             gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
            for (m=0; m < 3; m++) {
              for (n=0; n < 3; n++) {
                iter_gridinfo_w[gidy_to].disp[m][n] = iter_gridinfo_w[gidy_from].disp[m][n];
              }
            }
          }
        }
      }
    }
//     if (boundary[3].type==2) {
//     }
    //if((boundary[3].type == 2) && (workers_mpi.firsty || workers_mpi.lasty)) {
//     if (boundary[3].type == 2) {
//       copy_to   = boundary[3].points[0];
//       copy_from = boundary[3].proxy[0];
//       if (workers_mpi.firsty) {
//         y_first = boundary[3].proxy[0] - 1;
// //         printf("y_first=%ld\n", y_first);
//         constraint_uedge(boundary, y_first, x_start, x_end, copy_to, copy_from, -1.0);
//       }
//       else if (workers_mpi.lasty) {
//         y_last = boundary[3].proxy[0] + 1;
// //         printf("y_last=%ld\n", y_last);
//         constraint_uedge(boundary, y_last, x_start, x_end, copy_to, copy_from, 1.0);
//       }
//     }
  }
  if (strcmp(field_type, "F") == 0) {
   if(((boundary[4].type == 1) && (workers_mpi.firsty || workers_mpi.lasty)) || ((boundary[3].type == 3) && (workers_mpi.firsty && workers_mpi.lasty))) {
     for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[4].proxy[j];
        copy_to   = boundary[4].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from  = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to    = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;

            if (DIMENSION==2) {
              if (boundary[4].type ==1) {
                y_last = boundary[4].proxy[0] + 1;
                gridy_lbm_last = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_last;

                if (workers_mpi.lasty) {
                  lbm_gridinfo_w[gridy_lbm_last].fdis[6] = lbm_gridinfo_w[gridy_lbm_last].fdis[4];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[3] = lbm_gridinfo_w[gridy_lbm_last].fdis[1];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[7] = lbm_gridinfo_w[gridy_lbm_last].fdis[5];

                  for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }
                if (workers_mpi.firsty) {
                  y_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_first;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[5] = lbm_gridinfo_w[gridy_lbm_first].fdis[7];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[1] = lbm_gridinfo_w[gridy_lbm_first].fdis[3];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[4] = lbm_gridinfo_w[gridy_lbm_first].fdis[6];

                  for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
              }
            }
            
            if (DIMENSION==3){
              if (boundary[4].type ==1) {
                if (workers_mpi.lasty) {
                  // col fwd = col back
                  y_last = boundary[4].proxy[0] + 1;
                  gridy_lbm_last = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_last;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[4]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[3];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[8]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[7];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[12] =  lbm_gridinfo_w[gridy_lbm_last].fdis[11];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[13] =  lbm_gridinfo_w[gridy_lbm_last].fdis[14];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[18] =  lbm_gridinfo_w[gridy_lbm_last].fdis[17];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[20] =  lbm_gridinfo_w[gridy_lbm_last].fdis[19];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[22] =  lbm_gridinfo_w[gridy_lbm_last].fdis[21];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[23] =  lbm_gridinfo_w[gridy_lbm_last].fdis[24];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[26] =  lbm_gridinfo_w[gridy_lbm_last].fdis[25];

                  for (k=0;k<SIZE_LBM_FIELDS -4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }
                if (workers_mpi.firsty) {
                  // col back = colfwd
                  y_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_first;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[3]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[4];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[7]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[8];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[11] =  lbm_gridinfo_w[gridy_lbm_first].fdis[12]; // was 12 = 12
                  lbm_gridinfo_w[gridy_lbm_first].fdis[14] =  lbm_gridinfo_w[gridy_lbm_first].fdis[13];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[17] =  lbm_gridinfo_w[gridy_lbm_first].fdis[18];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[19] =  lbm_gridinfo_w[gridy_lbm_first].fdis[20];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[21] =  lbm_gridinfo_w[gridy_lbm_first].fdis[22];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[24] =  lbm_gridinfo_w[gridy_lbm_first].fdis[23];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[25] =  lbm_gridinfo_w[gridy_lbm_first].fdis[26];

                  for (k=0;k<SIZE_LBM_FIELDS -4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
              }
            }
            
            if (boundary[4].type == 3) {
              for (k=0;k<SIZE_LBM_FIELDS -4 ; k++) {
                lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gidy_from].fdis[k];
              }
              lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gidy_from].rho;
              for (dim=0; dim < DIMENSION; dim++) {
                lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gidy_from].u[dim];
              }
            }
          }
        }
      }
    }
    if((boundary[4].type == 2) && (workers_mpi.firsty || workers_mpi.lasty)) {
      for (j=0; j < 3; j++) { //Loop over three-buffer points
        copy_from = boundary[4].proxy[j];
        copy_to   = boundary[4].points[j];
        for (x=x_start; x<=x_end; x++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from  = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_from;
            gidy_to    = x*workers_mpi.layer_size + z*workers_mpi.rows_y + copy_to;
            if (DIMENSION ==2) {
              if (boundary[4].type ==2) {
                if (workers_mpi.lasty) {
                   y_last = boundary[4].proxy[0] + 1;
                   gridy_lbm_last = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_last;

                   lbm_gridinfo_w[gridy_lbm_last].fdis[3] = lbm_gridinfo_w[gridy_lbm_last].fdis[1] + 2.0*lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/3.0;

                   lbm_gridinfo_w[gridy_lbm_last].fdis[6] = lbm_gridinfo_w[gridy_lbm_last].fdis[4] - 0.5*(lbm_gridinfo_w[gridy_lbm_last].fdis[0]-lbm_gridinfo_w[gridy_lbm_last].fdis[2]) + lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/6.0;

                   lbm_gridinfo_w[gridy_lbm_last].fdis[7] = lbm_gridinfo_w[gridy_lbm_last].fdis[5] + 0.5*(lbm_gridinfo_w[gridy_lbm_last].fdis[0]-lbm_gridinfo_w[gridy_lbm_last].fdis[2]) + lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/6.0;

                   for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                   }

                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }
                if (workers_mpi.firsty) {
                  y_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x*workers_mpi.layer_size + z*workers_mpi.rows_y + y_first;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[1] = lbm_gridinfo_w[gridy_lbm_first].fdis[3] + 2.0*lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/3.0;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[4] = lbm_gridinfo_w[gridy_lbm_first].fdis[6] - 0.5*(lbm_gridinfo_w[gridy_lbm_first].fdis[0]-lbm_gridinfo_w[gridy_lbm_first].fdis[2]) + lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/6.0;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[5] = lbm_gridinfo_w[gridy_lbm_first].fdis[7] + 0.5*(lbm_gridinfo_w[gridy_lbm_first].fdis[0]-lbm_gridinfo_w[gridy_lbm_first].fdis[2]) + lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/6.0;

                  for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }

                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
              }
            }

          }
        }
      }
    }
  }
}
void copyYZ(struct bc_scalars *boundary, struct fields* gridinfo_w, char *field_type) {
  long gidy_from, gidy_to, y, a, k;
  long copy_from, copy_to;
  long j;
  long x, z;
  int m;
  if (strcmp(field_type, "PHI") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[0].type ==1) && (workers_mpi.firstx || workers_mpi.lastx)) || ((boundary[0].type == 3) && (workers_mpi.firstx && workers_mpi.lastx))) {
        copy_from = boundary[0].proxy[j];
        copy_to   = boundary[0].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = copy_from*workers_mpi.layer_size + z*workers_mpi.rows_y  + y;
            gidy_to            = copy_to*workers_mpi.layer_size   + z*workers_mpi.rows_y  + y;
            for (a=0; a < NUMPHASES; a++) {
              gridinfo_w[gidy_to].phia[a]     = gridinfo_w[gidy_from].phia[a];
              gridinfo_w[gidy_to].deltaphi[a] = gridinfo_w[gidy_from].deltaphi[a];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "MU") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[1].type ==1) && (workers_mpi.firstx || workers_mpi.lastx)) || ((boundary[1].type == 3) && (workers_mpi.firstx && workers_mpi.lastx))) {
        copy_from = boundary[1].proxy[j];
        copy_to   = boundary[1].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = copy_from*workers_mpi.layer_size  +  z*workers_mpi.rows_y + y;
            gidy_to            = copy_to*workers_mpi.layer_size    +  z*workers_mpi.rows_y + y;
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              gridinfo_w[gidy_to].compi[k] = gridinfo_w[gidy_from].compi[k];
              gridinfo_w[gidy_to].composition[k] = gridinfo_w[gidy_from].composition[k];
            }
          }
        }
      }
    }
    if((boundary[1].type ==2) && (workers_mpi.firstx || workers_mpi.lastx)) {
      if (workers_mpi.firstx) {
        x_first            = boundary[1].proxy[0] - 1;
        x_grid             = x_first;
      }
      if (workers_mpi.lastx) {
        x_last             = boundary[1].proxy[0] + 1;
        x_grid             = x_last;
      }
      for (j=0; j < 3; j++) {
        copy_to   = boundary[1].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gridx              = x_grid*workers_mpi.layer_size     +  z*workers_mpi.rows_y + y;
            gidy_to            = copy_to*workers_mpi.layer_size    +  z*workers_mpi.rows_y + y;
            Mu(ceq[NUMPHASES-1][NUMPHASES-1], Teq, NUMPHASES-1, gridinfo_w[gridx].compi);
            for (k=0; k < NUMCOMPONENTS-1; k++) {
              gridinfo_w[gidy_to].compi[k]       = gridinfo_w[gridx].compi[k];
              gridinfo_w[gidy_to].composition[k] = ceq[NUMPHASES-1][NUMPHASES-1][k];
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "T") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[2].type == 1) && (workers_mpi.firstx || workers_mpi.lastx)) || ((boundary[2].type == 3) && (workers_mpi.firstx && workers_mpi.lastx))) {
        copy_from = boundary[2].proxy[j];
        copy_to   = boundary[2].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = copy_from*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            gidy_to            = copy_to*workers_mpi.layer_size   + z*workers_mpi.rows_y + y;
            gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
          }
        }
      }
    }
  }
  if (strcmp(field_type, "U") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[3].type == 1) && (workers_mpi.firstx || workers_mpi.lastx)) || ((boundary[3].type == 3) && (workers_mpi.firstx && workers_mpi.lastx))) {
        copy_from = boundary[3].proxy[j];
        copy_to   = boundary[3].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from          = copy_from*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            gidy_to            = copy_to*workers_mpi.layer_size   + z*workers_mpi.rows_y + y;
//             gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
            for (m=0; m < 3; m++) {
              for (n=0; n < 3; n++) {
                iter_gridinfo_w[gidy_to].disp[m][n] = iter_gridinfo_w[gidy_from].disp[m][n];
              }
            }
          }
        }
      }
    }
  }
  if (strcmp(field_type, "F") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[4].type == 1) && (workers_mpi.firstx || workers_mpi.lastx)) || ((boundary[4].type == 3) && (workers_mpi.firstx && workers_mpi.lastx))) {
        copy_from = boundary[4].proxy[j];
        copy_to   = boundary[4].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from  = copy_from*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            gidy_to    = copy_to*workers_mpi.layer_size   + z*workers_mpi.rows_y + y;

            if (boundary[4].type == 3) {
              for (k = 0; k < SIZE_LBM_FIELDS-4 ; k++) {
                lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gidy_from].fdis[k];
              }
              lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gidy_from].rho;

              for (dim=0; dim < 3; dim++) {
                lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gidy_from].u[dim];
              }
              
            }

            if (DIMENSION == 2) {
              if (boundary[4].type == 1) {
                if (workers_mpi.firstx) {
                  x_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x_first*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
                  lbm_gridinfo_w[gridy_lbm_first].fdis[4] = lbm_gridinfo_w[gridy_lbm_first].fdis[6];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[0] = lbm_gridinfo_w[gridy_lbm_first].fdis[2];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[7] = lbm_gridinfo_w[gridy_lbm_first].fdis[5];

                  for (k=0; k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
                if (workers_mpi.lastx) {
                  x_last = boundary[4].proxy[0] + 1;
                  gridy_lbm_last = x_last*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[5] = lbm_gridinfo_w[gridy_lbm_last].fdis[7];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[2] = lbm_gridinfo_w[gridy_lbm_last].fdis[0];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[6] = lbm_gridinfo_w[gridy_lbm_last].fdis[4];

                  for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }
              }
            }
            
            //bounce back along x
            if (DIMENSION == 3){
              // printf("[BBX]") ;
              // LBM3D
              if(boundary[4].type == 1){
                //printf("\nFrom Boundary type 1 in YX plane along X direction") ;
                if (workers_mpi.firstx) {
                  //printf("\tFIRST x") ;
                  x_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x_first*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

                  /** the bounce back on YZ plane is coded as such:
                   * for 2D bounce back first x
                   * LHS = stream row fwd velocities 6,5,2 ; RHS = stream row back velocities 0,5,7
                   * So in 3D
                   * LHS = row fwd velocities and RHS is their correspondingly opposite velocities
                   * => row back velocities = row fwd velocities
                   */
                  
                  lbm_gridinfo_w[gridy_lbm_first].fdis[1]   =  lbm_gridinfo_w[gridy_lbm_first].fdis[2];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[7]   =  lbm_gridinfo_w[gridy_lbm_first].fdis[8];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[9]   =  lbm_gridinfo_w[gridy_lbm_first].fdis[10];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[13]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[14];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[15]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[16];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[19]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[20];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[21]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[22];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[23]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[24];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[26]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[25];

                  for (k=0;k<SIZE_LBM_FIELDS-4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
                if (workers_mpi.lastx) {
                  //printf("\tLAST x") ;
                  x_last = boundary[4].proxy[0] + 1;
                  gridy_lbm_last = x_last*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[2]   =  lbm_gridinfo_w[gridy_lbm_last].fdis[1];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[8]   =  lbm_gridinfo_w[gridy_lbm_last].fdis[7];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[10]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[9];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[14]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[13];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[16]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[15];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[20]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[19];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[22]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[21];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[24]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[23];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[25]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[26];

                  for (k=0;k<SIZE_LBM_FIELDS-4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }                
              }
            }
          }
        }
      }
      if((boundary[4].type == 2) && (workers_mpi.firstx || workers_mpi.lastx)) {
        copy_from = boundary[4].proxy[j];
        copy_to   = boundary[4].points[j];
        for (y=0; y < workers_mpi.rows_y; y++) {
          for (z=0; z < workers_mpi.rows_z; z++) {
            gidy_from  = copy_from*workers_mpi.layer_size + z*workers_mpi.rows_y + y;
            gidy_to    = copy_to*workers_mpi.layer_size   + z*workers_mpi.rows_y + y;
            if (DIMENSION == 2) {
              if (boundary[4].type == 2) {
                if (workers_mpi.firstx) {
                  x_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x_first*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[0] = lbm_gridinfo_w[gridy_lbm_first].fdis[2] + 2.0*lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/3.0;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[4] = lbm_gridinfo_w[gridy_lbm_first].fdis[6] - 0.5*(lbm_gridinfo_w[gridy_lbm_first].fdis[1]-lbm_gridinfo_w[gridy_lbm_first].fdis[3]) + lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/6.0;

                  lbm_gridinfo_w[gridy_lbm_first].fdis[7] = lbm_gridinfo_w[gridy_lbm_first].fdis[5] + 0.5*(lbm_gridinfo_w[gridy_lbm_first].fdis[1]-lbm_gridinfo_w[gridy_lbm_first].fdis[3]) + lbm_gridinfo_w[gridy_lbm_first].rho*boundary[4].value[j]/6.0;

                  for (k=0; k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }
                }
                if (workers_mpi.lastx) {
                  x_last = boundary[4].proxy[0] + 1;
                  gridy_lbm_last = x_last*workers_mpi.layer_size + z*workers_mpi.rows_y + y;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[2] = lbm_gridinfo_w[gridy_lbm_last].fdis[0] + 2.0*lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/3.0;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[6] = lbm_gridinfo_w[gridy_lbm_last].fdis[4] - 0.5*(lbm_gridinfo_w[gridy_lbm_last].fdis[1]-lbm_gridinfo_w[gridy_lbm_last].fdis[3]) + lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/6.0;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[5] = lbm_gridinfo_w[gridy_lbm_last].fdis[7] + 0.5*(lbm_gridinfo_w[gridy_lbm_last].fdis[1]-lbm_gridinfo_w[gridy_lbm_last].fdis[3]) + lbm_gridinfo_w[gridy_lbm_last].rho*boundary[4].value[j]/6.0;

                  for (k=0;k<9; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
void copyXY(struct bc_scalars *boundary, long x_start, long x_end, struct fields* gridinfo_w, char *field_type) {
  long gidy_from, gidy_to, y, a, k;
  long copy_from, copy_to;
  long j;
  long x, z;
  int m;
  if (strcmp(field_type, "PHI") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if(((boundary[0].type ==1) && (workers_mpi.firstz || workers_mpi.lastz)) || ((boundary[0].type == 3) && (workers_mpi.firstz && workers_mpi.lastz))) {
          copy_from = boundary[0].proxy[j];
          copy_to   = boundary[0].points[j];
          for (x=x_start; x <= x_end; x++) {
            for (y=0; y < workers_mpi.rows_y; y++) {
              gidy_from          = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
              gidy_to            = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y   + y;
              for (a=0; a < NUMPHASES; a++) {
                gridinfo_w[gidy_to].phia[a]     = gridinfo_w[gidy_from].phia[a];
                gridinfo_w[gidy_to].deltaphi[a] = gridinfo_w[gidy_from].deltaphi[a];
              }
            }
          }
        }
      }
   }
   if (strcmp(field_type, "MU") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
        if(((boundary[1].type ==1) && (workers_mpi.firstz || workers_mpi.lastz)) || ((boundary[1].type == 3) && (workers_mpi.firstz && workers_mpi.lastz))) {
          copy_from = boundary[1].proxy[j];
          copy_to   = boundary[1].points[j];
          for (x=x_start; x <= x_end; x++) {
            for (y=0; y < workers_mpi.rows_y; y++) {
              gidy_from          = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
              gidy_to            = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y   + y;
              for(k=0; k<NUMCOMPONENTS-1; k++) {
                gridinfo_w[gidy_to].compi[k]       = gridinfo_w[gidy_from].compi[k];
                gridinfo_w[gidy_to].composition[k] = gridinfo_w[gidy_from].composition[k];
              }
            }
          }
        }
      }
    }
   if (strcmp(field_type, "T") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if(((boundary[2].type ==1) && (workers_mpi.firstz || workers_mpi.lastz)) || ((boundary[2].type == 3) && (workers_mpi.firstz && workers_mpi.lastz))) {
          copy_from = boundary[2].proxy[j];
          copy_to   = boundary[2].points[j];
          for (x=x_start; x <= x_end; x++) {
            for (y=0; y < workers_mpi.rows_y; y++) {
              gidy_from                      = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
              gidy_to                        = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y   + y;
              gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
            }
          }
        }
      }
   }
   if (strcmp(field_type, "U") == 0) {
    for (j=0; j < 3; j++) { //Loop over three-buffer points
      if(((boundary[3].type ==1) && (workers_mpi.firstz || workers_mpi.lastz)) || ((boundary[3].type == 3) && (workers_mpi.firstz && workers_mpi.lastz))) {
          copy_from = boundary[3].proxy[j];
          copy_to   = boundary[3].points[j];
          for (x=x_start; x <= x_end; x++) {
            for (y=0; y < workers_mpi.rows_y; y++) {
              gidy_from                      = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
              gidy_to                        = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y   + y;
//               gridinfo_w[gidy_to].temperature = gridinfo_w[gidy_from].temperature;
              for (m=0; m < 3; m++) {
                for (n=0; n < 3; n++) {
                  iter_gridinfo_w[gidy_to].disp[m][n] = iter_gridinfo_w[gidy_from].disp[m][n];
                }
              }
            }
          }
        }
      }
   }
  if (strcmp(field_type, "F") == 0) {
   for (j=0; j < 3; j++) { //Loop over three-buffer points
     if(((boundary[4].type == 1) && (workers_mpi.firstz || workers_mpi.lastz)) || ((boundary[4].type == 3) && (workers_mpi.firstz && workers_mpi.lastz))) {
        copy_from = boundary[4].proxy[j];
        copy_to   = boundary[4].points[j];
        for (x=x_start; x <= x_end; x++) {
          for (y=0; y < workers_mpi.rows_y; y++) {
            gidy_from = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
            gidy_to = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y   + y;
            
            //bounce back along z
            if (DIMENSION == 3){
              if(boundary[4].type == 1){
                if (workers_mpi.firstz) {
                  z_first = boundary[4].proxy[0] - 1;
                  gridy_lbm_first = x*workers_mpi.layer_size + z_first*workers_mpi.rows_y + y;
                  
                  lbm_gridinfo_w[gridy_lbm_first].fdis[5]   =  lbm_gridinfo_w[gridy_lbm_first].fdis[6];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[9]   =  lbm_gridinfo_w[gridy_lbm_first].fdis[10];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[11]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[12];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[16]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[15];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[18]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[17];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[19]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[20];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[22]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[21];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[23]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[24];
                  lbm_gridinfo_w[gridy_lbm_first].fdis[25]  =  lbm_gridinfo_w[gridy_lbm_first].fdis[26];

                  for (k=0;k<SIZE_LBM_FIELDS-4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_first].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gridy_lbm_first].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_first].u[dim];
                  }

                }
                if (workers_mpi.lastz) {
                  z_last = boundary[4].proxy[0] + 1;
                  gridy_lbm_last = x*workers_mpi.layer_size + z_last*workers_mpi.rows_y + y;

                  lbm_gridinfo_w[gridy_lbm_last].fdis[6]   =  lbm_gridinfo_w[gridy_lbm_last].fdis[5];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[10]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[9];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[12]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[11];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[15]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[16];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[17]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[18];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[20]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[19];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[21]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[22];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[24]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[23];
                  lbm_gridinfo_w[gridy_lbm_last].fdis[26]  =  lbm_gridinfo_w[gridy_lbm_last].fdis[25];

                  for (k=0;k<SIZE_LBM_FIELDS-4; k++) {
                    lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gridy_lbm_last].fdis[k];
                  }
                  lbm_gridinfo_w[gidy_to].rho     = lbm_gridinfo_w[gridy_lbm_last].rho;
                  for (dim=0; dim < DIMENSION; dim++) {
                    lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gridy_lbm_last].u[dim];
                  }
                }                
              }
            }
          
          
            if (boundary[4].type == 3) {
              for (k=0;k<SIZE_LBM_FIELDS -4 ; k++) {
                lbm_gridinfo_w[gidy_to].fdis[k] = lbm_gridinfo_w[gidy_from].fdis[k];
              }
              lbm_gridinfo_w[gidy_to].rho = lbm_gridinfo_w[gidy_from].rho;
              for (dim=0; dim < DIMENSION; dim++) {
                lbm_gridinfo_w[gidy_to].u[dim] = lbm_gridinfo_w[gidy_from].u[dim];
              }
            }
          
          }
        }
      }
      // if((boundary[4].type == 2) && (workers_mpi.firstx || workers_mpi.lastx)) {
      //   copy_from = boundary[4].proxy[j];
      //   copy_to   = boundary[4].points[j];
      //   for (x=x_start; x <= x_end; x++) {
      //     for (y=0; y < workers_mpi.rows_y; y++) {
      //       gidy_from = x*workers_mpi.layer_size + copy_from*workers_mpi.rows_y + y;
      //       gidy_to = x*workers_mpi.layer_size + copy_to*workers_mpi.rows_y + y;



      //     }
      //   }
      // }// boundaty type 2
    }
  }




}



void apply_boundary_conditions(long taskid) {
  int i, j, field_num;
   if(boundary_worker) {
    //PERIODIC is the default boundary condition
    if (!((workers_mpi.firstx ==1) && (workers_mpi.lastx ==1))) {
      mpiboundary_left_right(taskid);
    }
    if (!((workers_mpi.firsty ==1) && (workers_mpi.lasty ==1))) {
      mpiboundary_top_bottom(taskid);
    }
    if (DIMENSION == 3) {
//       printf("firstz=%d, lastz=%d, rank_z=%d, front_node=%d, back_node=%d\n", workers_mpi.firstz, workers_mpi.lastz, workers_mpi.rank_z, workers_mpi.front_node, workers_mpi.back_node);
      if (!((workers_mpi.firstz ==1) && (workers_mpi.lastz ==1))) {
        mpiboundary_front_back(taskid);
      }
    }
    //PERIODIC is the default boundary condition

    for (i=0; i<6; i++) {
      if ((i==0) || (i==1)) {
        if (workers_mpi.firstx) {
          copyYZ(boundary[0], gridinfo_w, "PHI");
          copyYZ(boundary[0], gridinfo_w, "MU");
          if (!ISOTHERMAL) {
            copyYZ(boundary[0], gridinfo_w, "T");
          }
        }
        if (workers_mpi.lastx) {
          copyYZ(boundary[1], gridinfo_w, "PHI");
          copyYZ(boundary[1], gridinfo_w, "MU");
          if (!ISOTHERMAL) {
            copyYZ(boundary[1], gridinfo_w, "T");
          }
        }
      }
      if ((i==2) || (i==3)) {
//         copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "PHI");
//         copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "MU");
//         if (!ISOTHERMAL) {
//           copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "T");
//         }
        if (workers_mpi.lasty) {
          copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
          copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
          if (!ISOTHERMAL) {
            copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
          }
        }
        if (workers_mpi.firsty) {
          copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
          copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
          if (!ISOTHERMAL) {
            copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
          }
        }
      }
      if (DIMENSION == 3) {
//         if ((i==4)||(i==5)) {
//           copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "PHI");
//           copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "T");
//           }
//         }
        if ((i==4)||(i==5)) {
          if (workers_mpi.lastz) {
            copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
            copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
            if (!ISOTHERMAL) {
              copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
            }
          }
          if (workers_mpi.firstz) {
            copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
            copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
            if (!ISOTHERMAL) {
              copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
            }
          }
        }
      }
    }
  }
}
void apply_boundary_conditions_stress(long taskid) {
  int i, j, field_num;
   if(boundary_worker) {
    //PERIODIC is the default boundary condition
    if (!((workers_mpi.firstx ==1) && (workers_mpi.lastx ==1))) {
      mpiboundary_left_right_stress(taskid);
    }
    if (!((workers_mpi.firsty ==1) && (workers_mpi.lasty ==1))) {
      mpiboundary_top_bottom_stress(taskid);
    }
    if (DIMENSION == 3) {
//       printf("firstz=%d, lastz=%d, rank_z=%d, front_node=%d, back_node=%d\n", workers_mpi.firstz, workers_mpi.lastz, workers_mpi.rank_z, workers_mpi.front_node, workers_mpi.back_node);
      if (!((workers_mpi.firstz ==1) && (workers_mpi.lastz ==1))) {
        mpiboundary_front_back_stress(taskid);
      }
    }
    //PERIODIC is the default boundary condition

    for (i=0; i<6; i++) {
      if ((i==0) || (i==1)) {
        if (workers_mpi.firstx) {
//           copyYZ(boundary[0], gridinfo_w, "PHI");
//           copyYZ(boundary[0], gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyYZ(boundary[0], gridinfo_w, "T");
//           }
          copyYZ(boundary[0], gridinfo_w, "U");
        }
        if (workers_mpi.lastx) {
//           copyYZ(boundary[1], gridinfo_w, "PHI");
//           copyYZ(boundary[1], gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyYZ(boundary[1], gridinfo_w, "T");
//           }
          copyYZ(boundary[1], gridinfo_w, "U");
        }
      }
      if ((i==2) || (i==3)) {
//         copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "PHI");
//         copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "MU");
//         if (!ISOTHERMAL) {
//           copyXZ(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "T");
//         }
        if (workers_mpi.lasty) {
//           copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
//           copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
//           }
          copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "U");
        }
        if (workers_mpi.firsty) {
//           copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
//           copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
//           }
          copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "U");
        }
      }
      if (DIMENSION == 3) {
//         if ((i==4)||(i==5)) {
//           copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "PHI");
//           copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "MU");
//           if (!ISOTHERMAL) {
//             copyXY(boundary[i], workers_mpi.start[X], workers_mpi.end[X], gridinfo_w, "T");
//           }
//         }
        if ((i==4)||(i==5)) {
          if (workers_mpi.lastz) {
//             copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
//             copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
//             if (!ISOTHERMAL) {
//               copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
//             }
            copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "U");
          }
          if (workers_mpi.firstz) {
//             copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "PHI");
//             copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "MU");
//             if (!ISOTHERMAL) {
//               copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "T");
//             }
            copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "U");
          }
        }
      }
    }
  }
}

void apply_boundary_conditions_lbm(long taskid)
{
  int i, j, field_num;

  if(boundary_worker) {
    if (!( (workers_mpi.firstx ==1) && (workers_mpi.lastx ==1)) ) {
      mpiboundary_left_right_lbm(taskid);
    }
    if (!( (workers_mpi.firsty ==1) && (workers_mpi.lasty ==1)) ) {
      mpiboundary_top_bottom_lbm(taskid); 
    }
    if (DIMENSION == 3) {
      if (!( (workers_mpi.firstz ==1) && (workers_mpi.lastz ==1) )) {
        mpiboundary_front_back_lbm(taskid);
      }
    }
    

    for (i=0; i<6; i++) {
      if ((i==0) || (i==1)) {
        if (workers_mpi.firstx) {
          copyYZ(boundary[0], gridinfo_w, "F");
        }
        if (workers_mpi.lastx) {
          copyYZ(boundary[1], gridinfo_w, "F");
        }
      }
      if ((i==2) || (i==3)) {
        if (workers_mpi.lasty) {
          copyXZ(boundary[2], 0, workers_mpi.rows_x-1, gridinfo_w, "F");
        }
        if (workers_mpi.firsty) {
          copyXZ(boundary[3], 0, workers_mpi.rows_x-1, gridinfo_w, "F");
        }
      }
      if (DIMENSION == 3) {
        if ((i==4)||(i==5)) {
          if (workers_mpi.lastz) {
            copyXY(boundary[4], 0, workers_mpi.rows_x-1, gridinfo_w, "F");
          }
          if (workers_mpi.firstz) {
            copyXY(boundary[5], 0, workers_mpi.rows_x-1, gridinfo_w, "F");
          }
        }
      }
    }
  }
}




// void apply_boundary_conditions(long taskid){
// #ifdef PERIODIC
//   if (!((workers_mpi.firstx ==1) && (workers_mpi.lastx ==1))) {
//     mpiboundary_left_right(taskid);
//   }
// #endif
// #ifndef PERIODIC
//   if (workers_mpi.firstx ==1) {
//     copyYZ(0,5,gridinfo_w);
//     copyYZ(1,4,gridinfo_w);
//     copyYZ(2,3,gridinfo_w);
//   }
//   if (workers_mpi.lastx ==1) {
//     copyYZ(end[X]+1,end[X],gridinfo_w);
//     copyYZ(end[X]+2,end[X]-1,gridinfo_w);
//     copyYZ(end[X]+3,end[X]-2,gridinfo_w);
//   }
// #endif
// #ifdef PERIODIC_Y
//   if (!((workers_mpi.firsty ==1) && (workers_mpi.lasty ==1))) {
//      mpiboundary_top_bottom(taskid);
//   }
// //   copyXZ(2,MESH_Y-4,start,end,gridinfo_w);
// //   copyXZ(1,MESH_Y-5,start,end,gridinfo_w);
// //   copyXZ(0,MESH_Y-6,start,end,gridinfo_w);
// //
// //   copyXZ(MESH_Y-3,3,start,end,gridinfo_w);
// //   copyXZ(MESH_Y-2,4,start,end,gridinfo_w);
// //   copyXZ(MESH_Y-1,5,start,end,gridinfo_w);
// //
// #endif
// #ifdef ISOLATE_Y
//  if (workers_mpi.firsty ==1) {
//     copyXZ(2, 3, start[X], end[X], gridinfo_w);
//     copyXZ(1, 4, start[X], end[X], gridinfo_w);
//     copyXZ(0, 5, start[X], end[X], gridinfo_w);
//   }
//   if(workers_mpi.lasty ==1) {
//     copyXZ(rows_y-3, rows_y-4, start[X], end[X], gridinfo_w);
//     copyXZ(rows_y-2, rows_y-5, start[X], end[X], gridinfo_w);
//     copyXZ(rows_y-1, rows_y-6, start[X], end[X], gridinfo_w);
//   }
// #endif
// }
#endif
