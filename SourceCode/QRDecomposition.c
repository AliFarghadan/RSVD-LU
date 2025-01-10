
#include <slepcbv.h>
#include <Variables.h>

PetscErrorCode QRDecomposition(RSVD_vars *RSVD, Mat V)
{
	/*
		Performs QR decomposition on a given matrix 
	*/  

	PetscErrorCode        ierr;
	Mat                   Q_temp;
	BV                    Q;
	PetscInt              hh, mm, ss;
	PetscLogDouble        t1, t2;

	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** QR decomposition begins! ***\n");CHKERRQ(ierr);	

	ierr = BVCreateFromMat(V,&Q);CHKERRQ(ierr);
	ierr = BVSetType(Q,BVVECS);CHKERRQ(ierr);
	ierr = BVOrthogonalize(Q,NULL);CHKERRQ(ierr);
	ierr = BVCreateMat(Q,&Q_temp);CHKERRQ(ierr);
	ierr = MatCopy(Q_temp,V,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = BVDestroy(&Q);CHKERRQ(ierr);
	ierr = MatDestroy(&Q_temp);CHKERRQ(ierr);
	
	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** QR decomposition elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);	

	PetscFunctionReturn(0);

}



