
#include <petscmat.h>
#include <Variables.h>
#include <SetupFreqGrid.h>
#include <ReadUserInput.h>
#include <CreateRandomMat.h>
#include <CreateResultsDir.h>
#include <ReadWeightMats.h>
#include <SaveInputVarsCopy.h>

PetscErrorCode PreProcessing(RSVD_matrices *RSVDM, RSVD_vars *RSVD, Weight_matrices *Weight, Directories *dirs)
{

	/*
		Loading the operator, weight matrices, and creating the output directory
	*/  

	PetscErrorCode        ierr;
	PetscLogDouble        t1, t2;
	PetscInt              hh, mm, ss;
	PetscViewer           fd;

	PetscFunctionBeginUser;

	/*
		Reads in user input parameters
	*/

	ierr = ReadUserInput(RSVD, Weight, dirs);CHKERRQ(ierr);

	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*******************************************\n"
			"************** Problem info ***************\n*******************************************\n\n");CHKERRQ(ierr);

	/*
		Creates a folder for the results
	*/

	ierr = CreateResultsDir(dirs, "RSVDLU_ResolventModes_");CHKERRQ(ierr); 

	/*
		Saves a copy of the input variables in the results folder
	*/

	ierr = SaveInputVarsCopy(dirs);CHKERRQ(ierr);

	/*
		Initializes the time-stepping variables
	*/

	ierr = SetupFreqGrid(RSVD);CHKERRQ(ierr);

	/*
		Loads the operator
	*/

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Loading the operator: %s%s\n",dirs->RootDir,dirs->OperatorDir);CHKERRQ(ierr);

	ierr = PetscSNPrintf((char*)&dirs->IO_dir,PETSC_MAX_PATH_LEN,"%s%s",dirs->RootDir,dirs->OperatorDir);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,dirs->IO_dir,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&RSVDM->A_org);CHKERRQ(ierr);
	ierr = MatSetType(RSVDM->A_org,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatLoad(RSVDM->A_org,fd);CHKERRQ(ierr);
	ierr = MatGetSize(RSVDM->A_org,&RSVD->N,NULL);CHKERRQ(ierr);

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	if (RSVD->Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"Loading the operator elapsed time (N = %d) = %02d:%02d:%02d\n", (int)RSVD->N, (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	/*
		Reads weight and spatial matrices (if applicable)
	*/

	ierr = ReadWeightMats(RSVD, Weight, dirs);CHKERRQ(ierr);

	PetscFunctionReturn(0);
}


