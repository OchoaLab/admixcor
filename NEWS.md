# admixcor 0.0.0.9000 (2023-03-17)

First public version, includes:

* Exported functions `admixcor`, `align_Q`, and `rmsd_Q`
* Minimal unit tests

Missing documentation broadly.

# admixcor 0.0.1.9000 (2023-08-10)

- Function `admixcor` added option `L_algorithm` that allows potentially better optimizations of the internal `L` matrix (square root of Psi) using four algorithms:
  - `original` (default): performed unconstrained ridge regression on full L matrix (including lower triangle), then to solution set lower triangle to zero and caps remaining values between zero and one.
  - the three other options are variants of linear regression applied to proper L matrix (excluding lower triangle)
    - `nnls`: performs non-negative least squares, then sets caps values greater than 1 (no regularization)
	- `bvls`: performs bounded variable least squares, all L values bounded between 0 and 1 (no regularization)
	- `glmnet`: performs bounded ridge regression, all L values bounded between 0 and 1
- Added dependencies `nnls`, `bvls`, and `glmnet`.

Non-code updates:

* Added more detailed unit tests for the `update_{Q,R,L}` functions
* Refactored recurrent tests using new internal `validate_{Q,R,L}` functions.
* Added this `NEWS.md` file

# admixcor 0.0.2.9000 (2023-08-17)

- Function `admixcor` simplified internal `update_L` for all non-default `L_algorithm` cases so variables are optimized in independent blocks.  Original set them up as a single large problem to solve, which resulted in more complicated code and greater memory usage than the updated form.
  - Did have to write internal `glmnet_one_bounded` to solve the trivial case of bounded `glmnet` for a single variable, which `glmnet` itself refuses to solve (it requires two or more unknowns).

# admixcor 0.0.3.9000 (2023-08-18)

- Function `admixcor` added option `Q_algorithm` using the same four algorithms from `L_algorithm`.

# admixcor 0.0.4.9000 (2023-08-21)

- Function `admixcor` option `Q_algorithm` added algorithm "quadprog" to solve optimizations more exactly with all constraints into account.
  - Added dependency `quadprog`.

# admixcor 0.0.5.9000 (2023-08-21)

- Function `admixcor` option `Q_algorithm` added algorithm "quadprog-compact", hoping it performs better on account of sparsity of constraints

# admixcor 0.0.6.9000 (2023-09-08)

- Function `admixcor` now returns `ThetaSR` (the square root of the coancestry matrix `Theta`), `L` (the Cholesky decomposition of `Psi` that was already returned), and `R` (the rotation matrix relating `ThetaSR` to the non-negative admixture model decomposition).  This was done to ease troubleshooting, returning previously internal values that are needed to calculate on the outside the internal "square root" or "linearized" objective function.
- Function `objective` is now exported!  (Also for troubleshooting)
- Names of internal functions `initialize`, `theta_square_root`, and now exported function `objective` are now lowercased, so the style is consistent with the rest of the functions.
- Added minimal `README.md` with description of project and installation instructions.

# admixcor 0.0.7.9000 (2024-06-07)

- Function `admixcor` added option `kronecker` (default `FALSE` was previous behavior).  If `TRUE`, it frames the update of the internal `L` matrix more precisely using a Kronecker product (the previous formulation applied an approximate trick by pre-solving for one of the matrices, which works when there are no constraints but may be causing problems otherwise).  Only applies when `L_algorithm` is not "original".

# admixcor 0.0.8.9000 (2024-06-07)

- Function `admixcor` removed option `kronecker`.  The new default behavior is the previous `kronecker = TRUE` case, as evaluations confirmed the theoretical expectation that the Kronecker version is better than the previous approximate formulation, though performance differences are usually small.

# admixcor 0.0.9.9000 (2024-06-11)

- Function `admixcor` simplified option `L_algorithm`, which now defaults to `glmnet` and no longer accepts values `original` (performed worst) or `nnls` (too similar to `bvls`, which remains).  The default values for `gamma` and `delta` are now zero (both used to be 0.01).

# admixcor 0.0.10.9000 (2024-06-23)

- Function `admixcor` added option `vertex_refine`, which when `TRUE` (not default) explicitly tests vertices (all K unadmixed cases) for every individual, keeping the best if it improves upon previous solution at every update-Q iteration.

# admixcor 0.0.11.9000 (2024-07-01)

- Function `admixcor` added hack for `delta = 0` and `Q_algorithm %in% c('quadprog' 'quadprog-compact')` only, makes sure Psi is positive definite by ensuring that the internal L has diagonal values above 1e-7.  Before, especially if `gamma >= 1e-4`, could encounter these fatal errors:
  - Error in `quadprog::solve.QP(D, d, C, c, meq)`: matrix D in quadratic function is not positive definite!
  - Error in `quadprog::solve.QP.compact(D, d, Cmat, Cind, c, meq)`: matrix D in quadratic function is not positive definite!
- Internal: rearranged tests to see how often vertex-refine actually improved solutions (in a single step, always; but after a full round of iterations not often, those were commented out).  Also added comments on somewhat common failures that I haven't fixed yet

# admixcor 0.0.12.9000 (2024-07-01)

- Function `admixcor` increased minimum diagonal L value from 1e-7 to 1e-6 (there were still many problems on cluster runs)

# admixcor 0.0.13.9000 (2024-07-02)

- Function `admixcor` increased minimum diagonal L value again, from 1e-6 to 1e-5 (there were still many problems on cluster runs)

# admixcor 0.0.14.9000 (2024-07-02)

- Function `admixcor` increased minimum diagonal L value again, from 1e-5 to 1e-4 (there were still many problems on cluster runs)

# admixcor 0.0.15.9000 (2024-07-02)

- Function `admixcor` increased minimum diagonal L value again, from 1e-4 to 1e-3 (there were still many problems on cluster runs)

# admixcor 0.0.16.9000 (2024-07-02)

- Function `admixcor` increased minimum diagonal L value again, from 1e-3 to 1e-2 (there were still many problems on cluster runs)

# admixcor 0.0.17.9000 (2024-07-10)

- Function `admixcor` removed `Q_algorithm` cases `nnls`, `bvls`, `glmnet`, and `quadprog-compact`; in other words, only `quadprog` and `original` remain.

# admixcor 0.0.18.9000 (2024-09-09)

- Added functions `admixcor2` and `objective2`, which are based on optimizing `Psi` directly from the original objective (the previous `admixcor` and `objective` were based on the same "linearized" or "square root" objective used to optimize `Q`)
- Function `admixcor` added option `fix_L` that may fix a bug regarding the orientation of the internal `L` matrix (it appeared to be inconsistently transposed in certain steps).

# admixcor 0.0.19.9000 (2024-09-19)

- Actually exported functions `admixcor2` and `objective2`

# admixcor 0.0.20.9000 (2024-09-26)

- Added internal function `positive_definite`, which corrects intermediate `Psi` matrices to be positive definite, ultimately so that `chol` does not fail in the `admixcor2` function.  Corrects fatal errors in `admixcor2` such as:
  - Error in `chol.default(Psi1)`: the leading minor of order 3 is not positive

# admixcor 0.0.21.9000 (2024-10-18)

- Added functions `stretch_Q` that implements admixture stretching (a theoretical way to stretch admixture proportions to span as much of the simplex as possible without changing objective values), and `stretch_Psi` that applies the matching transformation to the subpopulations coancestry matrix.
- Removed package `nnls` from "Imports", since it was no longer being used).

# admixcor 0.0.22.9000 (2024-11-14)

- Function `admixcor` added options `stretch` to allow stretching in every iteration using function `stretch_Q`, `tol_stretch` to set the largest negative values allowed (temporarily within an interation), and `ties_none` to decide how to handle vertices associated with more than one ancestry.
- Function `stretch_Q` also added option `ties_none` as described above.  Internally this function now allows unstretched ancestries if no suitable vertices were identified.  The updates fix a bug that in some cases had this function hang (essentially stuck in an infinite loop) trying to choose unique vertices when there is no reasonable choice.

# admixcor 0.0.23.9000 (2024-11-22)

- Function `admixcor2` added option `stretch`, `tol_stretch`, and `ties_none` (see previous entry for `admixcor`).
- Function `admixcor` removed option `fix_L`, which is now always `TRUE`.
- Both functions `admixcor` and `admixcor2`:
  - Removed option `Q_type`, which is now always equal to `"random"`
  - Removed option `Q_algorithm`, which is now always equal to `"quadprog"`
  - Removed option `vertex_refine`, which is now always `FALSE`

# admixcor 0.0.24.9000 (2024-11-24)

- Functions `admixcor` and `admixcor2` removed options `L_algorithm` and `Psi_algorithm`, respectively.  These algorithm options now default to `"bvls"` when there is no penalization (`gamma = 0`, or `alpha = 0`, respectively), and default to `"glmnet"` otherwise.  Since `"bvls"` required no penalization before, the only new constraint is that now it is not possible to run `"glmnet"` without penalization, which was ill-advised anyway as it could result in errors or warnings.

# admixcor 0.0.25.9000 (2024-11-24)

- Functions `admixcor` and `admixcor2`:
  - Option `L_type` removed possible values `"identity"` and `"uniform"`, which in unit tests were associated with optimization problems (by initializing to overly symmetric solutions that do not favor one permutation over the others).
  - Internally changed order of updates (R, L, Q) when `stretch = TRUE`, to match non-stretching cases, which also seems to behave better avoiding strange errors.
  - Internal function `update_Q` now handles `delta = 0` case better, practically avoiding errors entirely in all cases.

# admixcor 0.0.26.9000 (2024-12-17)

- Functions `admixcor` and `admixcor2` removed option `stretch`, which now defaults to `TRUE`

# admixcor 0.0.27.9000 (2025-06-13)

- Functions `admixcor` and `admixcor2` added columns `L_singular` and `Q_error_rate` to output log table.

# admixcor 0.0.28.9000 (2025-06-13)

- Functions `admixcor` and `admixcor2` added option `restart_singular` to re-draw solutions if singularity is encountered instead the default that tries to regularize bad solutions and results in bad performance.

# admixcor 0.0.29.9000 (2025-06-14)

- Functions `admixcor` and `admixcor2` removed option `restart_singular` (it now defaults to TRUE) and removed column `Q_error_rate` to output log table.

# admixcor 0.0.30.9000 (2025-06-15)

- Functions `admixcor`, `admixcor2`, `objective` removed option `delta` (it now defaults to zero), and also `beta` in `admixcor2` and `objective2`.  Functions `objective` and `objective2` now return a vector with one less value (length 3 instead of 4).
