# Calculation codes for excitonic environmental (dielectric) effects in carbon nanotubes 

This repository contains several computer codes, mostly written in Fortran, used to perform calculations of exciton environmental effects, particularly dielectric constant effects, in single-walled carbon nanotubes (SWNTs).

## Contents of this README file

- [Main functionality of the program](#main-functionality-of-the-program)
- [First thing first: Check you can do make!](#first-thing-first-check-you-can-do-make)
- [Exciton energy (for $E_{ii}$ database or Kataura plot)](#exciton-energy-for-e_ii-database-or-kataura-plot)
- [Optimized $\kappa$](#optimized-kappa)
- [Least-Squares Regression of $\kappa$](#least-squares-regression-of-kappa)
- [Energy shift ($\Delta E_{\mathrm{env}}$)](#energy-shift-delta-e_mathrmenv)

---


## Main functionality of the program

The main computational task is the calculation of **exciton energies** for **many different dielectric constants** $\kappa$. Our implementation, published in [Appl. Phys. Lett. 97, 091905 (2010)](https://doi.org/10.1063/1.3485293), is basically based on the work by Jiang *et al.* [Phys. Rev. B **75**, 035407 (2007)] and Sato *et al.* [Phys. Rev. B **76**, 195446 (2007)].


The outputs of the main program are collected into a single database file. That database is then used by other programs for constructing the dielectric constant model $\kappa$ and for calculating the environmental energy shift $\Delta E_{\mathrm{env}}$.

All necessary programs are located in this GitHub repository, hereafter will be referred to as the `ROOT/` directory. More detailed usage instructions for each program might be found in the individual Fortran file or **Makefile** or 00README (if any) inside each subdirectory of `ROOT/`.

## First thing first: Check you can do make!

Each folder in this repository uses a `Makefile` to compile the Fortran source files and to remove intermediate build files after compilation. One representative example is:

```make
envkata: $(MODS_EXKATAMPI) $(OBJS_EXKATAMPI)
        mpif90 $(F90FLAGS) $(MODS_EXKATAMPI) $(OBJS_EXKATAMPI) \
        arguments-kappa.f90 envkata.f90 $(LAPACK) -o envkata.out

clean:
        rm -f *.a *.exe *.mod *.o *.out *~
````

### Target `envkata`

This target builds the executable for the program `envkata`. Its role is to compile and link all required module files, object files, and source files into the final executable:

* `envkata:`
  This is the name of the build target.

* `$(MODS_EXKATAMPI)`
  This variable usually contains the list of compiled module-related files or dependencies required before building `envkata`.

* `$(OBJS_EXKATAMPI)`
  This variable usually contains the object files needed for linking the final executable.

* `mpif90`
  This is the MPI Fortran compiler wrapper. It is used instead of plain `ifort` or `gfortran` because the code is intended to run in parallel using MPI.

* `$(F90FLAGS)`
  This variable stores compiler options, such as optimization flags, debugging flags, or standards settings.

* `arguments-kappa.f90 envkata.f90`
  These are the main Fortran source files explicitly compiled and linked in this target.
  In particular:

  * `arguments-kappa.f90` likely handles command-line or parameter input
  * `envkata.f90` contains the main program for exciton energy calculations

* `$(LAPACK)`
  This variable usually expands to the LAPACK and possibly BLAS libraries needed for numerical linear algebra routines.

* `-o envkata.out`
  This specifies the name of the generated executable. After successful compilation, the executable file will be `envkata.out`.

In practical terms, running

```bash
make envkata
```

will produce the executable `envkata.out`, provided all dependencies and compiler settings are correctly defined in the `Makefile`.

### Why `mpif90` is used

The program is designed for MPI-based parallel execution. Therefore, it must be compiled with an MPI-aware compiler wrapper. The generated executable can then be run, for example, with:

```bash
mpirun -np 4 ./envkata.out
```

where `4` is the number of MPI processes.

### Target `clean`

The `clean` target removes files generated during compilation:

```make
clean:
        rm -f *.a *.exe *.mod *.o *.out *~
```

Its purpose is to keep the repository directory tidy and to allow rebuilding from a clean state.

The removed file types are:

* `*.a`
  Static library files

* `*.exe`
  Executable files

* `*.mod`
  Fortran module files generated during compilation

* `*.o`
  Object files generated from source compilation

* `*.out`
  Output executables such as `envkata.out`

* `*~`
  Backup files, often created by text editors

The command

```bash
make clean
```

is useful when:

* recompiling everything from scratch
* avoiding inconsistencies from old object or module files
* cleaning the repository before sharing or archiving it

---

## Exciton energy (for $E_{ii}$ database or Kataura plot)

**Main program:** `ROOT/envkata/envkata.f90`  
**Database maker:** `ROOT/util/makeEii.f90`

Using `envkata.f90`, $E_{ii}$ energies are calculated for all $(n,m)$ SWNTs within

$$
0.5 < d_t < 3.0\ \mathrm{nm}
$$

and

$$
0 \leq \kappa \leq 8
$$

for $E_{11}$ up to $E_{44}$.

In a single run of the program, the required inputs are:

- a chirality value $(n,m)$,
- an index $i$ of $E_{ii}$,
- and a dielectric constant value $\kappa$.

Since a unified database is needed for many parameter combinations, the program is run repeatedly so that the desired diameter and $\kappa$ ranges are covered.

The diameter range is controlled by the variables `dtmin` and `dtmax` in `envkata.f90`, corresponding to the minimum and maximum diameter, respectively. In the present setup:

```fortran
dtmin = 0.5d0
dtmax = 3.0d0
```

The `envkata.f90` program requires MPI (Message Passing Interface) for parallel computation. It is therefore executed using batch scripts. An example batch script for calculating $E_{11}$ up to $E_{44}$ for semiconducting type-I SWNTs (`s1`) is shown below:

```bash
#!/bin/bash
for Eii in 1 2 3 4; do
    for x in 0 2 3 4 5 6 7; do
        for y in 0 1 2 3 4 5 6 7 8 9; do
            mpirun -np 5 ./envkata.out s1 $Eii $x$y
        done
    done
done
echo "finish all"
```

In the above script:

- `s1` denotes semiconducting type-I SWNTs,
- `xy` denotes $\kappa$.

With this notation:

- `xy = 10` corresponds to $\kappa = 1.0$,
- `xy = 20` corresponds to $\kappa = 2.0$,
- and so on.

To calculate $E_{ii}$ for metallic SWNTs and semiconducting type-II SWNTs, replace `s1` by `s0` and `s2`, respectively.

The `makeEii.f90` program is used to collect all separate data files generated by `envkata.f90` into a single database file named `eii.dat`.

The $E_{ii}$ energies and exciton size $l_k$ are arranged in arrays. Diameters and chiral angles are also stored in this database. The two main arrays are:

```fortran
Eii(n,m,i,nkappa,flagEii)
dttheta(n,m,flagdt)
```

where `nkappa` is an integer defined by

$$
\texttt{nkappa} = 10(\kappa - 1) + 1.
$$

In the `Eii` array:

- `flagEii = 1` returns $E_{ii}$,
- `flagEii = 2` returns the exciton size $l_k$.

In the `dttheta` array:

- `flagdt = 1` returns the diameter $d_t$,
- `flagdt = 2` returns the chiral angle $\theta$.

For example, to obtain $E_{11}$ for a $(6,5)$ SWNT with $\kappa = 2.2$, use:

```fortran
Eii(6,5,1,13,1)
```

The `Eii` and `dttheta` arrays are extensively used by the other programs. This avoids rerunning `envkata.f90` and therefore saves substantial computational time.

For a given $\kappa$, the `envkata.f90` program also outputs other constituents of $E_{ii}$:

- single-particle energy $\varepsilon$,
- exciton binding energy $E_{\rm bd}$,
- self-energy $\Sigma$,
- many-body energy $E_{\rm mb} = \Sigma - E_{\rm bd}$,

so that

$$
E_{ii} = \varepsilon + E_{\rm mb}.
$$

An example output is shown below for a single dielectric constant, $\kappa = 2.2$, and subband index $i = 1$ (that is, $E_{11}$) for semiconducting type-I SWNTs. This is only a small portion of the larger database.

| $n$ | $m$ | $d_t$ (nm) | $\theta$ (rad) | $E_{11}$ (eV) | $\epsilon$ (eV) | $\Sigma$ (eV) | $E_{\rm bd}$ (eV) | $E_{\rm mb}$ (eV) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 1 | 0.5252 | 0.1325 | 1.8666 | 2.7476 | 1.1223 | 0.8810 | 0.2413 |
| 6 | 4 | 0.6893 | 0.4086 | 1.3792 | 2.0554 | 0.9281 | 0.6762 | 0.2518 |
| 7 | 2 | 0.6492 | 0.2132 | 1.4886 | 2.2139 | 0.9602 | 0.7253 | 0.2349 |
| 7 | 5 | 0.8226 | 0.4277 | 1.1762 | 1.7379 | 0.7986 | 0.5618 | 0.2368 |
| 8 | 0 | 0.6356 | 0.0000 | 1.5544 | 2.2928 | 0.9665 | 0.7385 | 0.2281 |
| 8 | 3 | 0.7773 | 0.2669 | 1.2556 | 1.8559 | 0.8314 | 0.6003 | 0.2310 |
| 8 | 6 | 0.9564 | 0.4413 | 1.0250 | 1.5078 | 0.7033 | 0.4827 | 0.2205 |
| $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ |

*Example of the output format from the exciton energy program. Only the case of $\kappa = 2.2$ is shown here.*

## Optimized $\kappa$

**Main program:** `ROOT/diel/calkapp.f90`

The optimized $\kappa$ values are obtained by matching $E_{ii}$ for a given $(n,m)$ SWNT and subband index $i$ from experiments with the calculated $E_{ii}$ database.

Because it is generally impossible to satisfy exactly

$$
E_{ii}^{\rm exp} = E_{ii}^{\rm cal},
$$

we search for the optimized $\kappa$ that gives the smallest difference between experimental and calculated values.

The important parameters to save are:

- $n$,
- $m$,
- $p$,
- $i$,
- $E_{ii}$,
- the corresponding optimized $\kappa$.

This program uses the `Eii` and `dttheta` arrays generated previously by `makeEii.f90`.

Since the database uses a $\kappa$ step of only $0.1$, interpolation is used between neighboring $\kappa$ points in order to improve the optimized value. As a result, the output optimized $\kappa$ can be reported to four decimal places.

The input to the program is a set of:

- $(n,m)$ values,
- experimental $E_{ii}$ values,
- $p$ (cutting line index),
- $i$ (subband index).

An example input for an alcohol-assisted CVD sample is:

```text
#Alcohol-assisted cvd
 10  7  2.050  3  1
 11  0  1.560  2  2
 20  0  1.924  5  4
 17  1  2.292  4  3
```

## Least-Squares Regression of $\kappa$

**Main program:** `ROOT/diel/linkapp.f90`

The optimized $\kappa$ values obtained in the previous step are modeled using the functional form

$$
\kappa = C_{\kappa} \left[p^a \left(\frac{1}{d_t}\right)^b \left(\frac{1}{l_k}\right)^c \right] + C_x.
$$

For a particular sample, the best fit values of $(a,b,c)$, $C_{\kappa}$, and $C_x$ are determined by the least-squares method.

The result is then coupled with other experimental samples to check whether the same $(a,b,c)$ values can be used across samples. The best exponents with some available samples in our work in 2010, were found to be

$$
(a,b,c) = (0.80 \pm 0.10,\ 1.60 \pm 0.10,\ 0.40 \pm 0.05).
$$

The value of $C_{\kappa}$ for each sample is then recalculated by fixing a common $C_x$ value for all samples. This common $C_x$ determines the crossing point among the different fitting lines.

In this calculation, two conditions should be satisfied:

- the difference between experimental and calculated $E_{ii}$ obtained through $\kappa$ should be minimized,
- the linear regression correlation coefficient $R^2$ should be maximized.

With this treatment, the $\kappa$ values for $(E_{11}^S, E_{22}^S, E_{11}^M, E_{33}^S, E_{44}^S)$ in every experimental sample can be fitted by a single $\kappa$ function line, with only the $C_{\kappa}$ value changing from sample to sample.

Using `linkapp.f90`, the value of $C_{\kappa}$ can then be obtained from an input file of the same form as used in the optimized $\kappa$ calculation.

The normalized $C_{\kappa}$, denoted $\tilde{C}_{\kappa}$, is obtained by dividing the $C_{\kappa}$ of a given sample by that of the super-growth (SG) sample. With more experimental data available in the literature recently, we realized that the normalization of $C_{\kappa}$ might be performed using other SWNT samples, not necessarily the SG sample.


## Energy shift ($\Delta E_{\mathrm{env}}$)

**Main program:** `ROOT/diel/eshift.f90`

The ``environmental'' energy shift is expressed as

$$
\Delta E_{\mathrm{env}} = E_{ii}^{\rm SG} - E_{ii}^{\rm env}
\equiv \tilde{C}_{\kappa}
\left[A + B \left(\frac{p}{d_t}\right) + C \left(\frac{p}{d_t}\right)^2 \right].
$$

The coefficients $A$, $B$, and $C$ are obtained using multiple linear regression.

The above equation can be rearranged as

$$
\frac{\Delta E_{\mathrm{env}}}{\tilde{C}_{\kappa}} =
A + B \left(\frac{p}{d_t}\right) + C \left(\frac{p}{d_t}\right)^2,
$$

or written in the compact form

$$
Y = \beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3,
$$

where

- $X_1 = 1$,
- $X_2 = p/d_t$,
- $X_3 = (p/d_t)^2$

are the three independent variables, and

- $\beta_1 = A$,
- $\beta_2 = B$,
- $\beta_3 = C$

are the regression coefficients.

In matrix notation, if a particular sample contains $n$ data points for $\Delta E_{\mathrm{env}}$, then

$$
\vec{Y} = \vec{X}\,\vec{\beta},
$$

with

$$
\vec{Y} =
\begin{bmatrix}
y_1 \\
y_2 \\
\vdots \\
y_n
\end{bmatrix},
\qquad
\vec{X} =
\begin{bmatrix}
x_{11} & x_{12} & x_{13} \\
x_{21} & x_{22} & x_{23} \\
\vdots & \vdots & \vdots \\
x_{n1} & x_{n2} & x_{n3}
\end{bmatrix},
\qquad
\hat{\beta} =
\begin{bmatrix}
\beta_1 \\
\beta_2 \\
\beta_3
\end{bmatrix}.
$$

The coefficients $A$, $B$, and $C$ are then determined from

$$
\hat{\beta} = (\vec{X}^{T}\vec{X})^{-1}\vec{X}^{T}\vec{Y}.
$$