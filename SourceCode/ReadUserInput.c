
#include <petscmat.h>
#include <Variables.h>

PetscErrorCode ReadUserInput(RSVD_vars *RSVD, Weight_matrices *Weight, Directories *dirs)
{
	/*
		Reads in user input parameters
	*/

	PetscErrorCode        ierr;
	PetscBool             flg_set;
	char                  filename[PETSC_MAX_PATH_LEN];

	PetscFunctionBeginUser;

	ierr = PetscOptionsGetString(NULL, NULL,"-inputs",(char*)&filename,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify inputs variables via -inputs");CHKERRQ(ierr);
	ierr = PetscOptionsInsertFileYAML(PETSC_COMM_WORLD, NULL, filename, PETSC_FALSE);CHKERRQ(ierr);	
	ierr = PetscOptionsGetInt(NULL,NULL,"-RandSeed",&RSVD->RandSeed,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVD->RandSeed = 1373;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'RandSeed' variable not found. Setting 'RandSeed' to default value: %d\n", (int) RSVD->RandSeed);
	}	
	ierr = PetscOptionsGetInt(NULL,NULL,"-k",&RSVD->k,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'k'");CHKERRQ(ierr);
	} else if (RSVD->k < 1) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'k' must be a positive integer, current value: %d", (int) RSVD->k);CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetInt(NULL,NULL,"-q",&RSVD->q,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'q'");CHKERRQ(ierr);
	} else if (RSVD->q < 0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'q' must be a non-negative integer, current value: %d", (int) RSVD->q);CHKERRQ(ierr);
	}	
	ierr = PetscOptionsGetInt(NULL,NULL,"-Display",&RSVD->Display,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVD->Display = 2;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Display' variable not found. Setting 'Display' to default value: %d\n", (int) RSVD->Display);
	} else if (RSVD->Display > 2 || RSVD->Display < 0) {
		RSVD->Display = 2;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'Display' must be 0, 1 or 2. Setting 'Display' to default value: %d\n", (int) RSVD->Display);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-w_min",&RSVD->w_min,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'w_min'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-w_max",&RSVD->w_max,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'w_max'");CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(NULL,NULL,"-dw",&RSVD->dw,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'dw'");CHKERRQ(ierr);
	} else if (RSVD->dw <= 0) {
		SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'dw' must be positive, current value: %g", RSVD->dw);CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-DiscFlg",&RSVD->Disc.DiscFlg,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVD->Disc.DiscFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'DiscFlg' variable not found. Setting 'DiscFlg' to default value: %d\n", (int) RSVD->Disc.DiscFlg);
	}
	ierr = PetscOptionsGetReal(NULL,NULL,"-beta",&RSVD->Disc.beta,&flg_set);CHKERRQ(ierr);
	if (!flg_set && RSVD->Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Discounting flag is on! Either set 'DiscFlg' to zero or specify 'beta'");
	if (RSVD->Disc.beta < 0 && RSVD->Disc.DiscFlg) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"'beta' must be positive, current value: %g", RSVD->Disc.beta);
	ierr = PetscOptionsGetBool(NULL,NULL,"-TwoPI",&RSVD->TwoPI,&flg_set);CHKERRQ(ierr);
	if (!flg_set) {
		RSVD->TwoPI = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'TwoPI' variable not found. Setting 'TwoPI' to default value: %d\n", (int) RSVD->TwoPI);
	}
	ierr = PetscOptionsGetString(NULL, NULL,"-RootDir",(char*)&dirs->RootDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'RootDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-ResultsDir",(char*)&dirs->ResultsDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'ResultsDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL, NULL,"-OperatorDir",(char*)&dirs->OperatorDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
	if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OperatorDir'");CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,NULL,"-InvInputWeightFlg",&Weight->InvInputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight->InvInputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvInputWeightFlg' variable not found. Setting 'InputWeightFlg' to default value: %d\n", (int) Weight->InvInputWeightFlg);
	} else if (Weight->InvInputWeightFlg) {
		ierr = PetscOptionsGetString(NULL,NULL,"-InvInputWeightDir",(char*)&dirs->InvInputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvInputWeightDir'");CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-OutputWeightFlg",&Weight->OutputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight->OutputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputWeightFlg' variable not found. Setting 'OutputWeightFlg' to default value: %d\n", (int) Weight->OutputWeightFlg);
	} else if (Weight->OutputWeightFlg) {
		ierr = PetscOptionsGetString(NULL,NULL,"-OutputWeightDir",(char*)&dirs->OutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputWeightDir'");CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-InvOutputWeightFlg",&Weight->InvOutputWeightFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight->InvOutputWeightFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InvOutputWeightFlg' variable not found. Setting 'InvOutputWeightFlg' to default value: %d\n", (int) Weight->InvOutputWeightFlg);
	} else if (Weight->InvOutputWeightFlg) {
		ierr = PetscOptionsGetString(NULL,NULL,"-InvOutputWeightDir",(char*)&dirs->InvOutputWeightDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InvOutputWeightDir'");CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-InputMatrixFlg",&Weight->InputMatrixFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight->InputMatrixFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'InputMatrixFlg' variable not found. Setting 'InputMatrixFlg' to default value: %d\n", (int) Weight->InputMatrixFlg);
	} else if (Weight->InputMatrixFlg) {
		ierr = PetscOptionsGetString(NULL,NULL,"-InputMatrixDir",(char*)&dirs->InputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'InputMatrixDir'");CHKERRQ(ierr);
	}
	ierr = PetscOptionsGetBool(NULL,NULL,"-OutputMatrixFlg",&Weight->OutputMatrixFlg,&flg_set);CHKERRQ(ierr);	
	if (!flg_set) {
		Weight->OutputMatrixFlg = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: 'OutputMatrixFlg' variable not found. Setting 'OutputMatrixFlg' to default value: %d\n", (int) Weight->OutputMatrixFlg);
	} else if (Weight->OutputMatrixFlg) {
		ierr = PetscOptionsGetString(NULL,NULL,"-OutputMatrixDir",(char*)&dirs->OutputMatrixDir,PETSC_MAX_PATH_LEN,&flg_set);CHKERRQ(ierr);
		if (!flg_set) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Must specify 'OutputMatrixDir'");CHKERRQ(ierr);
	}
	
	PetscFunctionReturn(0);
}

