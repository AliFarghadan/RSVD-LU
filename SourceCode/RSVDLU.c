
#include <petscksp.h>
#include <Variables.h>
#include <CreateRandomMat.h>
#include <DirectAction.h>
#include <AdjointAction.h>
#include <PowerIteration.h>
#include <SVD4Response.h>
#include <SVD4Forcing.h>

PetscErrorCode RSVDLU(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight, Resolvent_matrices *Res, Directories *dirs, PetscInt iw)
{

	/*
		Performs the RSVD-LU algorithm to compute resolvent modes and gains for each freqency of interest
	*/  

	PetscErrorCode        ierr;
	KSP                   ksp;
	PC                    pc;
	Mat                   A;
	PetscInt              hh, mm, ss;
	PetscReal             w;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;
	
	/*
		Current frequency
	*/

	w = RSVD->w_min + iw * RSVD->dw;

	if (RSVD->Display) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"*******************************************\n"
			"       Processing iw = %d out of %d\n", (int)iw + 1, (int)RSVD->Nw);CHKERRQ(ierr);
		if (RSVD->TwoPI) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"    Computing resolvent modes at w = %g x 2pi\n", w/2/PETSC_PI);CHKERRQ(ierr);
		} else {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"    Computing resolvent modes at w = %g\n", w);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"*******************************************\n");CHKERRQ(ierr);
	}	

	/*
		Builds the reolvent operator
	*/

	ierr = MatDuplicate(RSVDM->A_org,MAT_COPY_VALUES,&A);CHKERRQ(ierr);
	ierr = MatScale(A, -1.);CHKERRQ(ierr);
	ierr = MatShift(A, PETSC_i * w);CHKERRQ(ierr);

	/*
		Discounting
	*/

	if (RSVD->Disc.DiscFlg) {
		ierr = MatShift(A,-RSVD->Disc.beta);CHKERRQ(ierr);
		if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---- Discounting with beta = %g ----\n\n", RSVD->Disc.beta);CHKERRQ(ierr);
	}

	/*
		Creates the initial random matrix of size Nxk
	*/

	ierr = CreateRandomMat(RSVDM, RSVD, dirs);CHKERRQ(ierr);

	/*
		Creates the KSP solver for R and R^*
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);

	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
	ierr = PCSetType(pc, PCLU);CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	ierr = KSPSetUp(ksp);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** LU decomposition elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*************************************************************************
		****************     RSVD - LU algorithm       *******************
		****************    for resolvent analysis     *******************
	**************************************************************************/

	/*
		Direct action
	*/

	ierr = DirectAction(ksp, RSVDM, RSVD, Weight);CHKERRQ(ierr);

	/*
		Power itertion
	*/	

	ierr = PowerIteration(ksp, RSVDM, RSVD, Weight);CHKERRQ(ierr);

	/*
		Reduced SVD to obtain response modes
	*/

	ierr = SVD4Response(RSVDM, RSVD, Weight, Res, dirs, iw);CHKERRQ(ierr);

	/*
		Adjoint action
	*/

	ierr = AdjointAction(ksp, RSVDM, RSVD, Weight);CHKERRQ(ierr);

	/*
		Reduced SVD to obtain forcing modes and gains
	*/	

	ierr = SVD4Forcing(RSVDM, RSVD, Weight, Res, dirs, iw);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Total elapsed time = %02d:%02d:%02d ***\n\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*
		Destroys variables
	*/

	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}	

