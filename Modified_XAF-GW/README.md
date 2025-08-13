# Computing partial chi using BerkeleyGW

**Reference:** J. Phys. Chem. Lett. 2021, 12, 36, 8841–8846

To compute partial chi, we need to generate a truncated wavefunction that 
contains only the relevant bands that are required to compute chi. 

This is achieved by asking `pw2bgw.x` to output the wavefunction only for a
specified set of band indices. There are two ways this code can perform
this operation:

1. Specifying an energy window around the Fermi level

   Flags:

   ```
   use_window = .true. (default: .false.)
   window = 0.5d0 (Width of energy window around Fermi level.
                   Units in eV)
   ```

2. Specifying the band indices to output the wavefunction

   Flags:

   ```
   use_bands = .true. (default: .false.)
   ct_band_min = 44 (min and max band indices to output 
   ct_band_max = 48  wavefunction)
   ```

`use_bands` and `use_window` cannot be both true at the same time.

Once the wavefunction is written out, you may check the correctness of
the output by using, for instance, `wfn_rho_vxc_info.x` in BerkeleyGW.

Then, the smaller wavefunction is used to compute chimat in `epsilon.x`.
chimat can then be added to other chimats using the rest of the XAF-GW
suite.

Currently, I provide the patch file for QE 6.4.1; modifications to more 
recent versions should be quite straightforward.

A simple example for monolayer graphene is included to illustrate the usage of
this code!

For more information on this method, please see the following reference and 
cite it if the method has been used for your work:

J. Phys. Chem. Lett. 2021, 12, 36, 8841–8846

Work done for this code was performed in the National University of Singapore (NUS).

Author: Nicholas Lin Quan Cheng (nlqcheng@nus.edu.sg)

Contact: Su Ying Quek (phyqsy@nus.edu.sg)
