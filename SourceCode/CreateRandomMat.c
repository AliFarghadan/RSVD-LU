
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode CreateRandomMat(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Directories *dirs)
{
	/*
		Generates a random matrix of size N x k
	*/  

	PetscErrorCode        ierr=0;
	PetscRandom           r;

	PetscFunctionBeginUser;

	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Generating a random forcing matrix\n\n");CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVDM->Y_hat);CHKERRQ(ierr);
	ierr = MatSetType(RSVDM->Y_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(RSVDM->Y_hat,PETSC_DECIDE,PETSC_DECIDE,RSVD->Nb,RSVD->k);CHKERRQ(ierr);
	ierr = MatSetUp(RSVDM->Y_hat);CHKERRQ(ierr);
	ierr = PetscRandomCreate(PETSC_COMM_WORLD,&r);CHKERRQ(ierr);
	ierr = PetscRandomSetSeed(r, RSVD->RandSeed);CHKERRQ(ierr);
	ierr = PetscRandomSeed(r);CHKERRQ(ierr);
	ierr = MatSetRandom(RSVDM->Y_hat,r);CHKERRQ(ierr);

	PetscFunctionReturn(0);
	
}

