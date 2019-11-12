# XAF-GW
The XAF-GW method enables GW calculations of large interface systems which can be separated into individual components without the formation of dangling bonds. However, the XAF-GW method can apply to layered systems with hybridization between the components, such as in the case of black phosphorus. 

XAF-GW stands for:

eXpand-chi
Add-chi
Full-sigma

The Add-chi method computes the chi matrix for the heterostructure (HS) as the sum of the chi matrices from individual components. This is shown to work even when bonding and anti-bonding states form between the separate components (JCTC).

The eXpand-chi method obtains the chi matrix of the supercell from the chi matrix of a unit cell.

The code provided here is based on BerkeleyGW/1.2.0, and routines for reading the chi matrices have been obtained from the BerkeleyGW code. However, this code is general and can be easily modified for other GW codes.

The BerkeleyGW/1.2.0 code is licensed under licensed under a free, open source, and permissive 3-clause modified BSD license.
BerkeleyGW, Copyright (c) 2011, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

Work done for this code was performed in the National University of Singapore (NUS). We ask that users of the XAF-GW method or of this code should cite the following reference:

J. Chem. Theory Comput. 2019, 15, 6, 3824-3835

If the BerkeleyGW code is used in conjunction with the XAF-GW method, please also cite the appropriate references for the BerkeleyGW code.

Author: Xuan Fengyuan (c2dxf@nus.edu.sg)

Contact: Su Ying Quek (phyqsy@nus.edu.sg) 
