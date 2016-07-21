##!/bin/bash
rm fort.*
#rm -rf tempdir_hdf5
#mpirun /home/nvarini/espresso_hdf5/espresso/PW/src/pw.x < si.scf_hdf5.in 2>&1|tee si.scf_hdf5.out
mpirun /home/nvarini/espresso_hdf5/espresso/PHonon/PH/ph.x < si.phX_hdf5.in 2>&1|tee si.phX_hdf5.out

#rm -rf tempdir_trunk
#mpirun /home/nvarini/espresso_hdf5/espresso_trunk/PW/src/pw.x < si.scf_trunk.in 2>&1|tee si.scf_trunk.out
#mpirun /home/nvarini/espresso_hdf5/espresso_trunk/PHonon/PH/ph.x < si.phX_trunk.in 2>&1|tee si.phX_trunk.out


#rm -rf tempdir_hdf5
#mpirun /home/nvarini/espresso_hdf5/espresso/PW/src/pw.x < si.scf_hdf5.in 2>&1|tee si.scf_hdf5.out
#mpirun /home/nvarini/espresso_hdf5/espresso/PHonon/PH/ph.x < si.phG_hdf5.in 2>&1|tee si.phG_hdf5.out

#rm -rf tempdir_trunk
#mpirun /home/nvarini/espresso_hdf5/espresso_trunk/PW/src/pw.x < si.scf_trunk.in 2>&1|tee si.scf_trunk.out
#mpirun /home/nvarini/espresso_hdf5/espresso_trunk/PHonon/PH/ph.x < si.phG_trunk.in 2>&1|tee si.phG_trunk.out

#rm -rf tempdir_hdf5
#mpirun /home/nvarini/espresso_hdf5/espresso/PW/src/pw.x < ch4.scf_hdf5.in 2>&1|tee ch4.scf.out
#mpirun /home/nvarini/espresso_hdf5/espresso/PHonon/PH/ph.x < ch4.nm_hdf5.in 2>&1|tee ch4.nm.out

