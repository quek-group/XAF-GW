Script `clean` wipes out all symbolic links and output files

1. Run `script_0` to create symbolic links

   **Note:** You may generate the k-grids using the script `generate_k_grids`
         if you're interested to know this. There is no need to run this
         as the k-grids are already in the input files.

2. Run `script_1` to get DFT wavefunctions and eigenvalues for the unshifted
   k-grid and shifted k-grid. Here, only the DFT wavefunctions that have energies
   within 0.5 eV of the Fermi energy are included in the 

   **Note:** A script `upgrade_to_hdf5` is provided for your convenience
         to convert wavefunctions to hdf5 format. If using hdf5
         wavefunctions, rememeber to add a flag `use_wfn_hdf5` later in
         epsilon.inp and sigma.inp

3. Run `script_2` to run epsilon.x using the DFT wavefunction with smaller subset
   of included bands. BerkeleyGW requires no modifications and runs as per usual. 

In my paper J. Phys. Chem. Lett. 2021, 12, 36, 8841–8846 I generated chimat file 
by using `skip_epsilon` and `dont_use_hdf5`, then I used machinery of XAF-GW 
(https://github.com/quek-group/XAF-GW) to do manipulation steps to chimat.
