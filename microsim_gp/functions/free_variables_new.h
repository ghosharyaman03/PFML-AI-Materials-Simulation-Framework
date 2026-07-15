#ifndef MICROSIM_FREE_VARS_HEADER
#define MICROSIM_FREE_VARS_HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include "global_vars.h"
#include "utility_functions_new.h"

void ms_FreeAllGlobalVariables()
{
    /**
     * Allocated in
     * ms_GlobalInit_PhaseFieldMatrices
     */

	free(start);
	free(end);
	free(averow);
	free(rows);
	free(offset);
	free(extra);

    ms_Free_3D(Diffusivity);
    ms_Free_3D(ceq);
    ms_Free_3D(cfill);
    ms_Free_3D(c_guess);
    ms_Free_3D(ceq_coeffs);
    ms_Free_3D(slopes);
    ms_Free_3D(dcbdT);
    ms_Free_3D(A);

    ms_Free_2D(DELTA_T);
    ms_Free_2D(DELTA_C);
    ms_Free_2D(dcbdT_phase);
    ms_Free_2D(B);
    ms_Free_2D(Beq);
    ms_Free_2D(dBbdT);
    ms_Free_1D(C);

    ms_Free_3D(cmu);
    ms_Free_3D(muc);
    
    ms_Free_4D(Rotation_matrix);
    ms_Free_4D(Inv_Rotation_matrix);
    ms_Free_1D(Rotated_qab);
    
    free(eigen_strain_phase);
    free(stiffness_phase);
    free(stiffness_phase_n);
    free(stiffness_t_phase);
    
    for (i = 0; i < 6; i++)
    {
        free(boundary[i]);
    }

    free(filling_type_phase);

    ms_Free_2D(Gamma);
    ms_Free_2D(tau_ab);
    ms_Free_2D(dab);
    ms_Free_2D(fab);

    ms_Free_3D(Gamma_abc);

    /**
     * Allocated in
     * ms_ReadInputParameters
     */

    ms_ArrayString_Free(Components, NUMCOMPONENTS);
    ms_ArrayString_Free(Phases,     NUMPHASES);
    ms_ArrayString_Free(Phases_tdb, NUMPHASES);
    ms_ArrayString_Free(phase_map,  NUMPHASES);

    if (ELASTICITY)
    {
        free(ec);
        free(e2);
        free(e4);
    }

    if (LBM)
    {
        free(rho_LBM);
        free(beta_c);
    }

    if (EIDT_SWITCH)
    {
        free(EIDT.comp_far_field);
        free(EIDT.comp_rate);
        free(EIDT.comp_center);
    }

    /**
     * Allocated in
     * initialize_variables
     */

    ms_Free_2D(dcdmu);
    ms_Free_3D(dcdmu_phase);
    ms_Free_2D(Ddcdmu);
    ms_Free_2D(inv_dcdmu);

    ms_Free_1D(deltamu);
    ms_Free_1D(deltac);
    ms_Free_1D(sum);
    ms_Free_1D(divphi);
    ms_Free_1D(lambda_phi);
    ms_Free_1D(divflux);
    ms_Free_1D(c_old);
    ms_Free_1D(c_new);
    ms_Free_1D(c_tdt);
    ms_Free_1D(divjat);

    /**
     * Allocated elsewhere
     */

    if (thermo_phase)
    {
        free(thermo_phase);
    }
}


#ifdef __cplusplus
}
#endif

#endif