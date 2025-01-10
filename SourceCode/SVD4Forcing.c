
#include <slepcsvd.h>
#include <Variables.h>

PetscErrorCode SVD4Forcing(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight, Resolvent_matrices *Res, Directories *dirs, PetscInt iw)
{
	/*
		Performs the economy SVD of matrices of size N \times k for Nw frequencies
	*/
	
	PetscErrorCode        ierr;
	PetscInt              ik, hh, mm, ss;
	PetscReal             sigma;
	Vec                   V;
	SVD                   svd;
	PetscViewer           fd;
	PetscLogDouble        t1, t2;


	PetscFunctionBeginUser;

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Reduced SVD begins! ***\n");CHKERRQ(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD,&Res->V_hat);CHKERRQ(ierr);
	ierr = MatSetType(Res->V_hat,MATDENSE);CHKERRQ(ierr);
	ierr = MatSetSizes(Res->V_hat,PETSC_DECIDE,PETSC_DECIDE,RSVD->N,RSVD->k);CHKERRQ(ierr);
	ierr = MatSetUp(Res->V_hat);CHKERRQ(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&Res->S_hat);CHKERRQ(ierr);
	ierr = VecSetSizes(Res->S_hat,PETSC_DECIDE,RSVD->k);CHKERRQ(ierr);
	ierr = VecSetUp(Res->S_hat);CHKERRQ(ierr);

	ierr = MatConjugate(RSVDM->Y_hat);CHKERRQ(ierr);

	ierr = SVDCreate(PETSC_COMM_WORLD,&svd);CHKERRQ(ierr);
	ierr = SVDSetOperators(svd,RSVDM->Y_hat,NULL);CHKERRQ(ierr);
	ierr = SVDSetDimensions(svd,RSVD->k,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
	ierr = SVDSolve(svd);CHKERRQ(ierr);

	for (ik=0; ik<RSVD->k; ik++) {
		ierr = MatDenseGetColumnVecWrite(Res->V_hat,ik,&V);CHKERRQ(ierr);
		ierr = SVDGetSingularTriplet(svd,ik,&sigma,V,NULL);
		ierr = VecSetValue(Res->S_hat,ik,sigma,INSERT_VALUES);
		ierr = MatDenseRestoreColumnVecWrite(Res->V_hat,ik,&V);CHKERRQ(ierr);
	}
	ierr = VecAssemblyBegin(Res->S_hat);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Res->S_hat);CHKERRQ(ierr);
	ierr = MatAssemblyBegin(Res->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Res->V_hat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	if (Weight->InvInputWeightFlg)  ierr = MatMatMult(Weight->W_f_sqrt_inv,Res->V_hat,MAT_REUSE_MATRIX,PETSC_DEFAULT,&Res->V_hat);CHKERRQ(ierr);
	
	/*
		Prints out the elapsed time
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Reduced SVD elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	ierr = SVDDestroy(&svd);CHKERRQ(ierr);
	ierr = MatDestroy(&RSVDM->Y_hat);CHKERRQ(ierr);

	/*
		Saves the forcing modes and gains
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Saving the forcing modes and gains begins! ***\n");CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"S_hat_iw",(int) iw+1,"_allK");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = VecView(Res->S_hat,fd);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s%d%s",dirs->FolderDir,"V_hat_iw",(int) iw+1,"_allK");CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	ierr = MatView(Res->V_hat,fd);CHKERRQ(ierr);

	ierr = MatDestroy(&Res->V_hat);CHKERRQ(ierr);
	ierr = VecDestroy(&Res->S_hat);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Saving the forcing modes and gains elapsed time = %02d:%02d:%02d ***\n\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	PetscFunctionReturn(0);

}



