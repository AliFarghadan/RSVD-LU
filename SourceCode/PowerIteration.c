
#include <petscksp.h>
#include <Variables.h>
#include <QRDecomposition.h>
#include <AdjointAction.h>
#include <DirectAction.h>

PetscErrorCode PowerIteration(KSP ksp, RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight)
{
	/*
		Performs power iteration for q times 
	*/

	PetscErrorCode        ierr;
	PetscInt              iq;

	PetscFunctionBeginUser;

	for (iq=0; iq<RSVD->q; iq++) {

		if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Inside power iteration, %d/%d *******\n\n",(int)iq+1,(int)RSVD->q);CHKERRQ(ierr);

		ierr = QRDecomposition(RSVD, RSVDM->Y_hat);CHKERRQ(ierr);
		ierr = AdjointAction(ksp, RSVDM, RSVD, Weight);CHKERRQ(ierr);
		ierr = MatConjugate(RSVDM->Y_hat);CHKERRQ(ierr);
		if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
		ierr = QRDecomposition(RSVD, RSVDM->Y_hat);CHKERRQ(ierr);
		ierr = DirectAction(ksp, RSVDM, RSVD, Weight);CHKERRQ(ierr);
		
	}

	if (RSVD->q > 0 && RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n******** Power iteration DONE! *************\n");CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



