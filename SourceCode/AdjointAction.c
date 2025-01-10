
#include <petscksp.h>
#include <Variables.h>
#include <ApplyWeightMats.h>

PetscErrorCode AdjointAction(KSP ksp, RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight)
{
	/*
		Computes R' \times \hat{F}, where (.)' indicates complex conjugate transpose
		For a modified resolvent operator, it computes B' * W_f_sqrt_inv' * R' * W_q_sqrt' * C' \times \hat{F} 
		In the latter case, the weight and input/output matrices are given as inputs
	*/

	PetscErrorCode        ierr;
	PetscLogDouble        t1, t2, t3, t4;
	Vec                   x;
	PetscInt              j, hh, mm, ss;

	PetscFunctionBeginUser;

	ierr = ApplyWeightMats(RSVDM, RSVD, Weight, 0, 1);CHKERRQ(ierr);

	ierr = MatConjugate(RSVDM->Y_hat);CHKERRQ(ierr);

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	for (j=0; j<RSVD->k; j++) {
		ierr = PetscTime(&t3);CHKERRQ(ierr);
		ierr = MatDenseGetColumnVecWrite(RSVDM->Y_hat, j, &x);CHKERRQ(ierr);
		ierr = KSPSolveTranspose(ksp, x, x);CHKERRQ(ierr);
		ierr = MatDenseRestoreColumnVecWrite(RSVDM->Y_hat, j, &x);CHKERRQ(ierr);
		ierr = PetscTime(&t4);CHKERRQ(ierr);
		hh   = (t4-t3)/3600;
		mm   = (t4-t3-3600*hh)/60;
		ss   = t4-t3-3600*hh-mm*60;
		if (RSVD->Display == 2) ierr = PetscPrintf(PETSC_COMM_WORLD,"Solving LU-decomposed system (adjoint) for mode number %d elapsed time = %02d:%02d:%02d ***\n", (int) j+1, (int)hh, (int)mm, (int)ss);CHKERRQ(ierr); 
	}
	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Solving LU-decomposed system (adjoint) for all modes elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);
	ierr = ApplyWeightMats(RSVDM, RSVD, Weight, 0, 0);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}



