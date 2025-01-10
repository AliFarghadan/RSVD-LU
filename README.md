
# RSVD-LU for resolvent analysis 

This repository contains the implementation of the Randomized Singular Value Decomposition (RSVD) tailored for resolvent analysis, referred to as **RSVD-LU**. The code computes resolvent modes over a user-defined range of frequencies. The computation process is detailed in our [paper](https://arxiv.org/pdf/2309.04617.pdf), drawing inspiration from the methodology in [Randomized Resolvent Analysis](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.033902). 

The notation aligns with our paper. The key distinction between this implementation and the original methodology lies in the treatment of response modes: we perform the Singular Value Decomposition (SVD) to obtain the response modes *prior* to the final adjoint action, skipping the QR decomposition and its subsequent recovery. This package is designed for seamless usage with PETSc and SLEPc libraries installed on your system.

---

## What are we computing?

We are using the RSVD-LU algorithm to compute resolvent modes of the linearized Navier-Stokes (LNS) equations. The resolvent operator $R$ is defined as

$$
R = C(\text{i}\omega I - A)^{-1} B,
$$

where $A$ is the LNS operator. $B$ and $C$ are input and output matrices, respectively. $R$ maps the input forcing to the output response in the frequency domain and can be expressed in terms of its singular vectors and values as 

$$
R = U \Sigma V^*.
$$

In our weighted formulation, we define the weighted resolvent operator $\tilde{R}$ as

$$
\tilde{R} = W_q^{1/2} C (\text{i}\omega I - A)^{-1} B W_f^{-1/2},
$$

where the resolvent modes are computed as

$$
\tilde{R} = \tilde{U} \Sigma \tilde{V}^*,
$$

$$
U = W^{-1/2}_q \tilde{U},
$$

$$
V = W^{-1/2}_f \tilde{V}.
$$

By the end of the simulation, resolvent modes (*i.e.*, gains, forcing, and response) are computed across all frequencies. In case $B$, $C$, $W_q^{-1/2}$, $W_q^{1/2}$, or $W_f^{-1/2}$ are not defined, we assume they are identity matrices.

---

## Installation and running RSVD-LU

### Installation

The installation process follows the steps outlined in the RSVD-Delta-t README, with variations in source codes and input variables due to algorithmic differences. Refer to the [RSVD-Delta-t README](https://github.com/AliFarghadan/RSVD-Delta-t/tree/Resolvent-analysis/README.md) for installation instructions. 

### Additional dependency: MUMPS

A key difference in compiling PETSc/SLEPc for RSVD-LU is the inclusion of the **MUMPS** package, required for parallel computation of LU decomposition for sparse matrices. Below is a recommended PETSc configuration:

```bash
./configure --with-debugging=0 --with-scalar-type=complex
--with-64-bit-indices PETSC_ARCH=complex-opt --download-mumps
```

**Note**: You may need to modify `RSVDLU.c` function the source code if you have a more suitable solver for your problem.

## List of input variables

Here is a list of variables used for our the RSVD-LU algorithm. Please note that you will need to modify these variables to suit your own projects.

```yaml

# Root directory (string)
# This directory must exist before the simulation
RootDir:            /path/to/root/directory/ 

# Results directory (string)
# The resolvent modes/gains will be saved in RootDir/ResultsDir/RSVDLU_ResolventModes_<int>
# This directory must exist before the simulation
ResultsDir:         /path/to/results/ 

# The linearized operator (string)
# The operator directory is defined as RootDir/OperatorDir
# This directory must exist before the simulation, OperatorDir can be same as ResultsDir
OperatorDir:        /path/to/A # expected a matrix saved in binary format

# Number of test vectors (integer)
k:                  10

# Number of power iterations (integer)
q:                  1

# Minimum frequency to resolve (real)
w_min:              -1.00

# Maximum frequency to resolve (real)
w_min:              1.00

# frequency rsolution (real > 0)
dw:                 0.05

# Convert frequencies to angular frequencies (boolean)
# if true: w <-- 2*pi*w, otherwise, w <-- w
TwoPI:              false

# Display options (integer 0 <= Display <= 2)
Display:            2

# Discounting flag (boolean)
# Applies discounting for unstable linear systems, A <-- A - beta I
DiscFlg:            false

# beta value when "DiscFlg = True", otherwise beta is ignored (real > 0)
beta:               0.1

# Seeding random number to replicate data if needed (integer)
RandSeed:           14

# Inverse input weight flag (boolean)
# The inverse input weight directory must exist before the simulation if true, otherwise ignored
InvInputWeightFlg:  false
# Inverse input weight directory when "InvInputWeightFlg = True" (string)
# This directory is defined as RootDir/InvInputWeightDir
InvInputWeightDir:  /path/to/W_f_sqrt_inv # expected a matrix saved in binary format

# Output weight flag (boolean)
# The input weight directory must exist before the simulation if true, otherwise ignored
OutputWeightFlg:    false
# Output weight directory when when "OutputWeightFlg = True" (string)
# This directory is defined as RootDir/OutputWeightDir
OutputWeightDir:    /path/to/W_q_sqrt # expected a matrix saved in binary format

# Inverse output weight flag (boolean)
# The inverse input weight directory must exist before the simulation if true, otherwise ignored
InvOutputWeightFlg: false
# Inverse output weight directory when "InvOutputWeightFlg = True" (string)
# This directory is defined as RootDir/InvOutputWeightDir
InvOutputWeightDir: /path/to/W_q_sqrt_inv # expected a matrix saved in binary format

# Input matrix flag (boolean)
# The input matrix directory must exist before the simulation if true, otherwise ignored
InputMatrixFlg:     false
# Input matrix directory when when "InputMatrixFlg = True" (string)
# This directory is defined as RootDir/InputMatrixDir
InputMatrixDir:     /path/to/B # expected a matrix saved in binary format

# Output matrix flag (boolean)
# The output matrix directory must exist before the simulation if true, otherwise ignored
OutputMatrixFlg:    false
# Output matrix directory when when "OutputMatrixFlg = True" (string)
# This directory is defined as RootDir/OutputMatrixDir
OutputMatrixDir:    /path/to/C # expected a matrix saved in binary format

```

### Important notes:

- For boolean flags, the following values are equivalent:
  - `False`, `false`, and `0` all represent `false`
  - `True`, `true`, and `1` all represent `true`

- Mathematical expressions in inputs are not allowed; for example, `2 x 2`, `√3`, and `π` will generate errors.

- For integer values, ensure only integers are provided. Decimal values such as `2.5`, `1.`, or `3.14` will cause errors.

- When using discounting, ensure `beta` is positive.

- Error messages will be displayed, and the simulation will terminate if an error occurs.

### Prerequisite modules

Make sure you have the following prerequisite modules loaded/installed:

- **OpenMPI** (or another MPI package such as MPICH)
- A **C++ compiler**
- **PETSc** and **SLEPc** packages

---

### RSVD-LU variables

- `RootDir`: Specifies the root directory path for the simulation.
- `ResultsDir`: Defines the path to the results directory where output files will be saved. This directory must exist within `RootDir`. If it does not, the system will create the directory at the specified path. Ensure you have write access to the root directory, or an error message will be displayed.
- `OperatorDir`: Specifies the directory path for the linearized operator matrix. If the operator is located in the `RootDir`, you only need to provide the operator name (e.g., `A_GL`). Otherwise, specify the relative path to the operator from `RootDir` (e.g., `matrices/A_GL`).
- `k`: The number of test vectors.
- `q`: The number of power iterations.
- `w_min`: The minimum frequency to resolve.
- `w_max`: The maximum frequency.
- `dw`: Frequency step size.
- `TwoPI`: Indicates whether to convert frequencies to angular frequencies. If `true`, the output frequencies are converted to angular frequencies by multiplying by $2\pi$.
- `InputForcingDir`: Defines the path to the input forcing $(\hat{F})$ directory.
- `InvInputWeightFlg`: Determines whether the inverse input weight is used.
- `InvInputWeightDir`: Defines the path to the inverse input weight $(W_f^{-1/2})$ directory.
- `InputMatrixFlg`: Determines whether the input matrix is used.
- `InputMatrixDir`: Defines the path to the input matrix $(B)$ directory.
- `InvOutputWeightFlg`: Determines whether the inverse output weight is used.
- `InvOutputWeightDir`: Defines the path to the inverse output weight $(W_q^{-1/2})$ directory.
- `OutputWeightFlg`: Determines whether the output weight is used.
- `OutputWeightDir`: Defines the path to the output weight $(W_q^{1/2})$ directory.
- `OutputMatrixFlg`: Determines whether the output matrix is used.
- `OutputMatrixDir`: Defines the path to the output matrix $(C)$ directory.
- `Display`: Controls the amount of information printed during computation, ranging from 0 (no output) to 2 (verbose output):
  - `Display = 0`: Minimal output with no information displayed.
  - `Display = 1`: Standard output, displaying:
    - problem information
    - elapsed time of the LU decomposition at each frequency
    - Elapsed time of QR and SVDs
    - Elapsed time of saving modes
    - Estimated remaining time
  - `Display = 2`: Detailed output, including everything from `Display = 1`, plus the elapsed time of solving LU-decomposed system for every test vector.
- `RandSeed`: Indicates the seed number for random number generation. Using the same number of cores and `RandSeed` value allows for repeatable results in simulations.
- `DiscFlg`: A boolean flag indicating whether to use a discounting strategy.
- `beta`: Specifies the `beta` value when `DiscFlg = true`. It is ignored if `DiscFlg = false`.

---

### Saving resolvent modes and gains

- A folder is created in the results directory with a fixed prefix `RSVDLU_ResolventModes_<int>`, where `<int>` is an integer starting from 0. If `RSVDLU_ResolventModes_i` exists, the code increments the integer until a unique folder name is found, ensuring that results from different simulations are not overwritten.

- For each frequency, response modes (each of size `N × k`) are saved as `U_hat_iw<int>_allK`, where `<int>` represents the integer index of the frequency. Similarly, forcing modes (each of size `N × k`) are saved as `V_hat_iw<int>_allK`. The corresponding gains, containing `k` singular values for each frequency, are saved as `S_hat_iw<int>_allK` of size `k × 1`.
- The indices correspond to frequencies within $\Omega = \omega_{\text{min}}:dw:\omega_{\text{max}}$ starting with 1.
- For instance, `U_hat_iw1_allK`, `V_hat_iw1_allK`, and `S_hat_iw1_allK` contain the response, forcing, and gains, respectively, associated with the first frequency ($\omega$ = `w_min`).

## Practical recommendation

For real-valued matrices, the resolvent modes are symmetric around $\omega = 0$. Hence, you can set `w_min = 0` without losing generality.

### References

* [Scalable resolvent analysis for three-dimensional flows](https://arxiv.org/pdf/2309.04617.pdf), *Journal of Computational Physics*, 2024
* [Randomized resolvent analysis](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.033902), *Physical Review Fluids*, 2020
  
### Contact information

Ali Farghadan, University of Michigan, aliii@umich.edu\
Aaron Towne, University of Michigan, towne@umich.edu

### Cite as

```cite
@article{farghadan2024scalable,
  title={Scalable resolvent analysis for three-dimensional flows},
  author={Farghadan, Ali and Martini, Eduardo and Towne, Aaron},
  journal={Journal of Computational Physics},
  pages={113695},
  year={2024},
  publisher={Elsevier}
}
```

