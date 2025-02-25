# VASPBERRY
Berry curvature and Chern number calculations with the output (WAVECAR) of VASP code.
VASPBERRY is written for the post-processing purpose of the VASP outputs, i.e., WAVECAR the Bloch wavefunction information. VASPBERRY can compute Berry curvature and Chern number via Fukui's method [See J. Phys. Soc. Jap. 74, 1674 (2005)]. In addition Circular dichroism also can be evaluated. Since it directly reads the wavefunction coefficients, one can also obtain real space wavefunction character psi_nk(r) by simple command.

# Download Git version
* git clone --branch master  https://github.com/Infant83/TBFIT.git

# Compile
* ~~Serial version~~ : 
    > ~~ifort -fpp -assume byterecl -mkl -o vaspberry vaspberry.f~~
* Multicore version : 
    > mpiifort vaspberry.f90 -fpp -assume byterecl -I$MKLROOT/include/intel64/lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -o vaspberry

* ~~Note for gfortran~~:
    ~~For gfortran, please use vaspberry_gfortran_serial.f for the compilation. This only support non-parallel calculations.~~
    ~~For the compilation, for example~~
    > ~~gfortran -L/usr/local/lib/lapack/ -l lapack -o vaspberry vaspberry_gfortran_serial.f~~

# Features
* Berry curvature calculation
* Compute Chern number for certain band(s) which are well-isolated over the BZ
* Circular dichroism (optical selectivity response to the circulary polarized light)
* Wavefunction plot (Gamma point only in the current version)

# Usage
* Instruction and possible options
> ./vaspberry -h
* Berry curvature calculation and Chern number (ex, k-grid: 12x12, multiband berry curvature from 1-th to 18-th band)
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18
* Circular dichroism [ex, transition rate from 11-th to 12-th state by right(+) polarized light]
> ./vaspberry -kx 12 -ky 12 -cd 1 -ii 11 -if 12
* Real space wavefunction plot [ex, to plot 18-th state with 1-st k-point (if it is gamma point), with 40x40x40 grid for density file]
> ./vaspberry -wf 18 -k 1 -ng 40,40,40
> NONE: current version only support gamma point for wavefunction plot. (there is some problem in boundary region in off-gamma k-point)
* If your system is semimetallic, there can be following error messages: "error. !!! ne(k) /= ne(k') !!!". This is due to that the number of occupied states for certain k-point (ne(k)) counted based on the calculated Fermi level is differ over the Brillouin zone. In this case, one can explicitly specify the number of electrons (NE) of your system, so that VASPBERRY calculate berry curvature with "NE" bands. 
> ./vaspberry -kx 12 -ky 12 -ii 1 -if 18 -ne 18

# Example
* 1H-MoS2 : Berry curvature and Chern number
* Quantum Anomalous Hall effect (Trypheny-lead lattice) : See H.-J. Kim, C. Li, J. Feng, J.-H. Cho, and Z. Zhang, PRB 93, 041404(R) (2016) (the example files will be provided upon request)
* Circular dichroism : See S.-W. Kim, H.-J. Kim, S. Cheon, and T.-H. Kim, Phys. Rev. Lett. accepted (2021) (the example will be provided upon reasonable request).

# Contributors
* Hyun-Jung Kim: Main developer (h.kim@fz-juelich.de, PGI-1/IAS-1, Forschungszentrum Jülich)
* Sun-Woo Kim: Circular dichroism, Kubo formular (kimsunwoo821@gmail.com, Department of Physics, Sungkyunkwan University)

# Citation of the code:
@software{Kim_VASPBERRY_2018,author = {Kim, Hyun-Jung},doi = {10.5281/zenodo.1402593},month = {8},title = {{VASPBERRY}},url = {https://github.com/Infant83/VASPBERRY},version = {1.0},year = {2018}}

@article{PhysRevLett.128.046401,
  title = {Circular Dichroism of Emergent Chiral Stacking Orders in Quasi-One-Dimensional Charge Density Waves},
  author = {Kim, Sun-Woo and Kim, Hyun-Jung and Cheon, Sangmo and Kim, Tae-Hwan},
  journal = {Phys. Rev. Lett.},
  volume = {128},
  issue = {4},
  pages = {046401},
  numpages = {6},
  year = {2022},
  month = {Jan},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.128.046401},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.128.046401}
}


