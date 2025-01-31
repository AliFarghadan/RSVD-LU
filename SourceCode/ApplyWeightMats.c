
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ApplyWeightMats(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight, PetscBool DirAdj, PetscBool before)
{
	/*
		Applies weight, input and output matrices if defined
	*/

	PetscErrorCode        ierr=0;
	Mat                   Y;

	PetscFunctionBeginUser;

	if (DirAdj) { // direct 
		if (before) { // forcing 
			if (Weight->InvInputWeightFlg)  {
				ierr = MatMatMult(Weight->W_f_sqrt_inv,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVDM->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight->InputMatrixFlg)  {
				ierr = MatMatMult(Weight->B,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVD->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response 
			if (Weight->OutputMatrixFlg) {
				ierr = MatMatMult(Weight->C,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight->OutputWeightFlg) {
				ierr = MatMatMult(Weight->W_q_sqrt,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVDM->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}		
		}
	} else { // adjoint
		if (before) { // forcing
			if (Weight->OutputWeightFlg) {
				ierr = MatHermitianTranspose(Weight->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight->W_q_sqrt);CHKERRQ(ierr);
				ierr = MatMatMult(Weight->W_q_sqrt,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVDM->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight->W_q_sqrt, MAT_INPLACE_MATRIX, &Weight->W_q_sqrt);CHKERRQ(ierr);
			}
			if (Weight->OutputMatrixFlg) {
				ierr = MatHermitianTranspose(Weight->C, MAT_INPLACE_MATRIX, &Weight->C);CHKERRQ(ierr);
				ierr = MatMatMult(Weight->C,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight->C, MAT_INPLACE_MATRIX, &Weight->C);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
		} else { // response
			if (Weight->InputMatrixFlg)  {
				ierr = MatHermitianTranspose(Weight->B, MAT_INPLACE_MATRIX, &Weight->B);CHKERRQ(ierr);
				ierr = MatMatMult(Weight->B,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatHermitianTranspose(Weight->B, MAT_INPLACE_MATRIX, &Weight->B);CHKERRQ(ierr);
				ierr = MatDestroy(&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDuplicate(Y,MAT_COPY_VALUES,&RSVDM->Y_hat);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);
			}
			if (Weight->InvInputWeightFlg)  {
				ierr = MatHermitianTranspose(Weight->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight->W_f_sqrt_inv);CHKERRQ(ierr);
				ierr = MatMatMult(Weight->W_f_sqrt_inv,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Y);CHKERRQ(ierr);
				ierr = MatCopy(Y,RSVDM->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
				ierr = MatDestroy(&Y);CHKERRQ(ierr);	
				ierr = MatHermitianTranspose(Weight->W_f_sqrt_inv, MAT_INPLACE_MATRIX, &Weight->W_f_sqrt_inv);CHKERRQ(ierr);
			}
		}
	}

	PetscFunctionReturn(0);
	
}



