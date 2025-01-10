
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode SetupFreqGrid(RSVD_vars *RSVD)
{
	/*
		Initializes the number of frequencies to resolve based on user input
	*/

	PetscErrorCode        ierr=0;

	PetscFunctionBeginUser;

	if (RSVD->w_max < RSVD->w_min) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"w_min = %g must be smaller than w_max = %g", RSVD->w_min, RSVD->w_max);CHKERRQ(ierr); 
	RSVD->Nw = PetscFloorReal((RSVD->w_max - RSVD->w_min) / RSVD->dw) + 1;

	RSVD->w_min *= RSVD->TwoPI ? 2*PETSC_PI : 1;
	RSVD->w_max *= RSVD->TwoPI ? 2*PETSC_PI : 1;
	RSVD->dw    *= RSVD->TwoPI ? 2*PETSC_PI : 1;

	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,\
		"@ k          = %d\n@ q          = %d\n@ w_min      = %g\n@ w_max      = %g\n@ dw         = %g\n\n",\
		(int)RSVD->k,(int)RSVD->q,RSVD->w_min,RSVD->w_max,RSVD->dw);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}




