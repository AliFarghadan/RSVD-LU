
#ifndef VARIABLES_H
#define VARIABLES_H

typedef struct {
	PetscBool       DiscFlg;                                /* discounting flag */
	PetscReal       beta;                                   /* discounting parameter */
} Discounting;

typedef struct {
	PetscInt        N;                                      /* problem size (state dimension) */
	PetscInt        Nb;                                     /* input size */
	PetscInt        Nc;                                     /* output size */
	PetscInt        k;                                      /* number of test vectors */
	PetscInt        q;                                      /* number of power iterations */
	PetscInt        Nw;                                     /* number of input/output frequencies to resolve */
	PetscInt        RandSeed;                               /* seeding random number to replicate data if desired */
	PetscReal       w_min;                                  /* min frequency */
	PetscReal       w_max;                                  /* max frequency */
	PetscReal       dw;                                     /* frequency resolution */
	PetscBool       TwoPI;                                  /* base frequency multiplies by 2*pi if true */
	PetscBool       RealOperator;                           /* real-valued matrix if true, otherwise complex-valued */
	Discounting     Disc;                                   /* discounting variables */
	PetscInt        Display;                                /* 0: None, 1: Partial, 2: Full progress display */
} RSVD_vars;

typedef struct {
	Mat             W_q_sqrt;                               /* output weight matrix */
	Mat             W_q_sqrt_inv;                           /* inverse output weight matrix */
	Mat             W_f_sqrt_inv;                           /* inverse input weight matrix */
	Mat             B;                                      /* input matrix */
	Mat             C;                                      /* output matrix */
	PetscBool       InvInputWeightFlg;                      /* inverse input weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool       OutputWeightFlg;                        /* output weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool       InvOutputWeightFlg;                     /* inverse output weight matrix from the specified directory if true, otherwise identity matrix */
	PetscBool       InputMatrixFlg;                         /* input matrix from the specified directory if true, otherwise identity matrix */
	PetscBool       OutputMatrixFlg;                        /* output matrix from the specified directory if true, otherwise identity matrix */
} Weight_matrices;

typedef struct {
	Mat             A_org;                                  /* LNS operator */
	Mat             Y_hat;                                  /* RSVD matrix */
} RSVD_matrices;

typedef struct {
	Mat             U_hat;                                  /* response resolvent modes */
	Mat             V_hat;                                  /* forcing resolvent modes */
	Vec             S_hat;                                  /* resolvent gains */
} Resolvent_matrices;

typedef struct {
	char            RootDir[PETSC_MAX_PATH_LEN];            /* root directory */
	char            ResultsDir[PETSC_MAX_PATH_LEN];         /* results folder */
	char            OperatorDir[PETSC_MAX_PATH_LEN];        /* LNS operator directory */
	char            filename[PETSC_MAX_PATH_LEN];           /* filename */
	char            IO_dir[PETSC_MAX_PATH_LEN];             /* I/O directory */
	char            FolderDir[PETSC_MAX_PATH_LEN];          /* results folder directory */
	char            InvInputWeightDir[PETSC_MAX_PATH_LEN];  /* inverse input weight directory */ 
	char            OutputWeightDir[PETSC_MAX_PATH_LEN];    /* output weight directory */ 
	char            InvOutputWeightDir[PETSC_MAX_PATH_LEN]; /* inverse output weight directory */
	char            InputMatrixDir[PETSC_MAX_PATH_LEN];     /* input matrix directory */
	char            OutputMatrixDir[PETSC_MAX_PATH_LEN];    /* output matrix directory */
} Directories;

#endif

