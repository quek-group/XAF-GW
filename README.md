# XAF-GW
The high computational cost of GW calculation for large 2D interface system is partially due to the polarizability(chi) matrix calculation.
We first cut the total system into layers and then cut each layer into unit cell. In Reference: JCTC we demonstrate that the total 
chi matrix can be written as the sum of each layer's chi matrix. Each layer's chi in a super cell can be expanded from it unit cell.

This reporsitory contains the codes based on BerkeleyGW/1.2.0 to expand unit cell chi matrix to supercell 
and summation of the supercell chi.

For instructions see NOTE.

Author: Xuan Fengyuan (c2dxf@nus.edu.sg)
