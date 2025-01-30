
#include <slepcsvd.h>
#include <Variables.h>

PetscErrorCode SVD4Response(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight, Resolvent_matrices *Res, Directories *dirs, PetscInt iw)
{
	/*
		Performs the economy SVD of a matrix of size N x k
		We perform SVD instead of QR to obtain U 
		This is generally more accurate than performing QR and recovering it later
	*/
	
	PetscErrorCode        ierr;
	PetscInt              ik, hh, mm, ss;
	Mat                   Y;
	Vec                   U;
	SVD                   svd;
	PetscViewer           fd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Y);CHKERRQ(ierr);
	ierr = MatSetType(Y,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Y,PETSC_DECIDE,PETSC_DECIDE,RSVD->Nc,RSVD->k);CHKERRQ(ierr);
	ierr = MatSetUp(Y);CHKERRQ(ierr);

	ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
	ierr = SVDSetOperators(svd,RSVDM->Y_hat,NULL);CHKERRQ(ierr);
	ierr = SVDSetDimensions(svd,RSVD->k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = SVDSolve(svd);CHKERRQ(ierr);

	for (ik=0; ik<RSVD->k; ik++) {
		ierr = MatDenseGetColumnVecWrite(Y,ik,&U);CHKERRQ(ierr);
		ierr = SVDGetSingularTriplet(svd,ik,NULL,U,NULL);
		ierr = MatDenseRestoreColumnVecWrite(Y,ik,&U);CHKERRQ(ierr);
	}
	
	ierr = SVDDestroy(&svd);CHKERRQ(ierr);

	ierr = MatAssemblyBegin(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Y,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatCopy(Y,RSVDM->Y_hat,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
	ierr = MatDestroy(&Y);CHKERRQ(ierr);

	/*
		Prints out the elapsed time
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Reduced SVD elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*
		Saves the response modes
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Saving the response modes begins! ***\n");CHKERRQ(ierr);

	if (Weight->InvOutputWeightFlg) ierr = MatMatMult(Weight->W_q_sqrt_inv,RSVDM->Y_hat,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&Res->U_hat);CHKERRQ(ierr);
	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"U_hat_iw",(int) iw+1,"_allK");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res->U_hat,fd);CHKERRQ(ierr);

	ierr = MatDestroy(&Res->U_hat);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVDM->Y_hat);CHKERRQ(ierr);

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Saving the response modes elapsed time = %02d:%02d:%02d ***\n\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



