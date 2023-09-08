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
- Added minimal `README.md` with description of project and installation instructions.

