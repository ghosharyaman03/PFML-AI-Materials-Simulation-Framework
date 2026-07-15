#ifndef MICROSIM_UTILITY_2_HEADER
#define MICROSIM_UTILITY_2_HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "global_vars.h"
#include "utility_functions_new.h"

void IO_Initialize()
{
    IO_NUM_FIELDS = NUMPHASES;

    if ((FUNCTION_F!=5) && (!GRAIN_GROWTH)) {
        IO_NUM_FIELDS += 2*(NUMCOMPONENTS-1);
    }

    if (ELASTICITY) {
        if (DIMENSION == 2) {
            IO_NUM_FIELDS += 2;
        } else {
            IO_NUM_FIELDS += 3;
        }
    }

    if (LBM) {
        IO_NUM_FIELDS += 1;
        if (DIMENSION == 2) {
            IO_NUM_FIELDS += 2 + 9;
        } else {
            IO_NUM_FIELDS += 3 + 27;
        }
    }

    if(!ISOTHERMAL) {
        IO_NUM_FIELDS += 1;
    }

    IO_FIELD_NAMES = msMalloc(char*, IO_NUM_FIELDS);
    for (size_t i = 0; i < IO_NUM_FIELDS; i++) {
        IO_FIELD_NAMES[i] = msMalloc(char, 64);
    }

    size_t field_id = 0 ;
    size_t k;

    for (k = 0; k < NUMPHASES; k++) {
        sprintf(IO_FIELD_NAMES[field_id++], "/%s",Phases[k]);
    }

    if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
        for (k = 0; k < (NUMCOMPONENTS-1); k++) {
            sprintf(IO_FIELD_NAMES[field_id++],"/Mu_%s", Components[k]);
        }
        for (k = 0; k < (NUMCOMPONENTS-1); k++) {
            sprintf(IO_FIELD_NAMES[field_id++], "/Composition_%s", Components[k]);
        }
    }

    if (ELASTICITY) {
        sprintf(IO_FIELD_NAMES[field_id++],"/Ux");
        sprintf(IO_FIELD_NAMES[field_id++],"/Uy");

        if (DIMENSION == 3) {
            sprintf(IO_FIELD_NAMES[field_id++],"/Uz");
        }
    }

    if (LBM) {
        sprintf(IO_FIELD_NAMES[field_id++],"/LBM_Ux");
        sprintf(IO_FIELD_NAMES[field_id++],"/LBM_Uy");

        if (DIMENSION == 3) {
            sprintf(IO_FIELD_NAMES[field_id++],"/LBM_Uz");
        }
    }

    if (!ISOTHERMAL) {
        sprintf(IO_FIELD_NAMES[field_id++],"/T");
    }

    if(LBM) {
        sprintf(IO_FIELD_NAMES[field_id++],"/LBM_rho");
        for(k = 0 ; k < 9; k ++) {
            sprintf(IO_FIELD_NAMES[field_id++],"/LBM_f_%ld", k);   
        }
        if (DIMENSION == 3) {
            for(k = 9 ; k < 27; k ++)
            {
                sprintf(IO_FIELD_NAMES[field_id++],"/LBM_f_%ld", k);   
            }
        }
    }

    assert(field_id == IO_NUM_FIELDS && "VALUES MUST MATCH");

    if(taskid == MASTER)
    {
        fprintf(stderr, "\n- [%d] IO_NUM_FIELDS = %ld",taskid, IO_NUM_FIELDS) ;
        for(k = 0 ; k < IO_NUM_FIELDS; k++)
        {
            fprintf(stderr, "\n- [%d] IO_FIELD_NAMES[%ld] = `%s`",taskid, k, IO_FIELD_NAMES[k]) ;
        }
    }
}

void IO_MallocGlobal()
{
    long index_count = workers_mpi.rows[X]*workers_mpi.rows[Y];
    if (DIMENSION == 3)
    {
        index_count *= workers_mpi.rows[Z];
    }
    IO_BUFFER = ms_Malloc_2D(IO_NUM_FIELDS, index_count);
}

void IO_FreeGlobal()
{
    ms_Free_2D(IO_BUFFER);
    for (size_t i = 0; i < IO_NUM_FIELDS; i++)
    {
        free(IO_FIELD_NAMES[i]);
    }
    free(IO_FIELD_NAMES);
}

#ifdef __cplusplus
}
#endif

#endif
