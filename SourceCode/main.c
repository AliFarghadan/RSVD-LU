
/* 	
	Author: Ali Farghadan, 2025
	Produced at the Univeristy of Michigan, Towne Lab
	Reference papers: 
	1. Scalable resolvent analysis for three-dimensional flows, JCP, 2024
	2. Randomized resolvent analysis, PRF, 2020
*/

/*  RSVD-LU algorithm to compute resolvent modes for a range of frequencies */

/*

	List of inputs *** Description ************************************  Format

	RSVD variables:
	A                  linear (or LNS) operator                          matrix
	B                  input matrix as defined in the reference paper 1  matrix
	C                  output matrix as defined in the reference paper 1 matrix
	W_q_sqrt           W_q^(1/2) as defined in the reference paper 1     matrix
	W_q_sqrt_inv       W_q^(-1/2) as defined in the reference paper 1    matrix
	W_f_sqrt_inv       W_f^(-1/2) as defined in the reference paper 1    matrix
	k                  number of test vectors                            integer
	q                  number of power iterations                        integer
	w_min              min frequency                                     real 
	w_max              max frequency                                     real 
	dw                 frequency reolution                               real > 0
	TwoPI              base frequency multiplies by 2*pi if true         boolean
	RootDir            root directory                                    string
	ResultsDir         results directory (RootDir/ResultsDir)            string
	beta               beta value for discounting (A <-- A - beta I)     real > 0
	RandSeed           seeding random number                             integer
	DiscFlg            applies discounting for unstable linear systems   boolean
	InputMatrixFlg     applies input matrix                              boolean 
	OutputMatrixFlg    applies output matrix                             boolean 
	InputWeightFlg     applies input weight matrix                       boolean
	InvInputWeightFlg  applies inverse input weight matrix               boolean 
	InvOutputWeightFlg applies inverse output weight matrix              boolean 
	Display            display options                                   integer
	    case 1) Display = 0: nothing
	    case 2) Display = 1: problem information + elapsed time of the LU decomposition at each frequency
	    case 3) Display = 2: "Display = 1" information + elapsed time of solving LU-decomposed system for all test vectors
	SaveResultsOpt     saving resolvent modes options                   integer
	    case 1) SaveResultsOpt = 1: saves resolvent modes as k  matrices of size N x Nw
	    case 2) SaveResultsOpt = 2: saves resolvent modes as Nw matrices of size N x k

	List of outputs ** Description ************************************  Format

	U                  response resolvent modes                          matrix
	V                  forcing resolvent modes                           matrix
	Sigma              resolvent gains                                   vector

*/

/* 	
	List of input libraries and functions
*/

#include <slepcsys.h>
#include <Variables.h>
#include <PreProcessing.h>
#include <RSVDLU.h>

/* 	
	Beginning of the simulation
*/

int main(int argc,char **args)
{

	/* 	
		Defines variable types
	*/

	PetscErrorCode        ierr;                             /* Petsc error code */
	Directories           dirs;                             /* I/O directories */
	RSVD_vars             RSVD;                             /* RSVD variables */
	Weight_matrices       Weight;                           /* weight and input/output matrices */
	RSVD_matrices         RSVDM;                            /* RSVD matrices */
	Resolvent_matrices    Res;                              /* resolvent modes and gains */
	PetscLogDouble        t1, t2;                           /* timing variables */
	PetscInt              hh, mm, ss;                       /* timing variables */

	/*
		Initializes the SLEPc
	*/

	ierr = SlepcInitialize(&argc,&args,(char*)0,NULL);if (ierr) return ierr;

	/*
		Reads input vaiables, creates and reads in the required matrices before running the algorithm
	*/

	ierr = PreProcessing(&RSVDM, &RSVD, &Weight, &dirs);CHKERRQ(ierr);
	
	/*************************************************************************
		******************     RSVD - LU algorithm     *******************
		******************    for resolvent analysis   *******************
	**************************************************************************/

	if (RSVD.Display) ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*******************************************\n"
			"***************** RSVD-LU *****************\n*******************************************\n\n");CHKERRQ(ierr);

	ierr = PetscTime(&t1);CHKERRQ(ierr);
	
	for (PetscInt iw=0; iw<RSVD.Nw; iw++) {

		ierr = RSVDLU(&RSVDM, &RSVD, &Weight, &Res, &dirs, iw);CHKERRQ(ierr);

	}

	/*
		Prints out the elapsed time and exits
	*/

	ierr = PetscTime(&t2);CHKERRQ(ierr);
	hh   = (t2-t1)/3600;
	mm   = (t2-t1-3600*hh)/60;
	ss   = t2-t1-3600*hh-mm*60;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"The results directory: %s\n\n",dirs.FolderDir);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"DONE :))\n\n*** Entire simulation elapsed time = %02d:%02d:%02d ***\n", (int)hh, (int)mm, (int)ss);CHKERRQ(ierr);

	ierr = PetscOptionsClear(NULL);CHKERRQ(ierr);
	ierr = SlepcFinalize();
	return ierr;

}

/* 	
	The end!
*/



