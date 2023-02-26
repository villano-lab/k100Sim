# k100Sim -- Geant4 simulation for the K100 cryogenic setup at UMN

This code was previously used with Geant4 v10.1.2 but not documented very well. We have migrated
to v10.7.4 and ROOT 6.24.2; both substantially advanced versions from what we started with, but
still not cutting edge (G4 is up to v11.1.1 and ROOT is up to v6.28.0 as of this writing). 


#Installing with G4 v10.7.4, ROOT 6.24.2, Ubuntu 20, gcc 9.4.0, standard=cxx14.


##Geant4 Install:

 - Download the source code: [v10.7.4](https://geant4.web.cern.ch/download/10.7.4.html#releasenotes)
 - make a build directory.
 - inside the build directory do:

```
cmake -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_INSTALL_DATA=ON -DCMAKE_CXX_STANDARD=14 -DCMAKE_PREFIX_PATH="/home/villaa/install/geant4/geant4-v10.7.4-install/lib/cmake/Geant4/" -DCMAKE_INSTALL_PREFIX=/home/villaa/install/geant4/geant4-v10.7.4-install /home/villaa/install/geant4/geant4-v10.7.4
```

Note the flags, if you miss the visualization, it will not be able to do visualization, if you
miss the install data flag it won't install the data (this happens at the make install step). If
you miss the forcing of the `cxx14` standard it will compile with a standard that might be
inconsistent with ROOT. The compilation will fail, hard. 

 - also inside the build directory do `$make -jN`, N is the number of parallel processors 
 - then do `$make install` (check that data directories are being downloaded!!)

 - make sure after you're installed to set the environmental variables like:

`source /home/villaa/install/geant4/geant4-v10.7.4-install/bin/geant4.sh` <br>
`export G4WORKDIR=/home/userrname/geant4bin`

The last one is to ensure that the `make install` process puts the binaries in a central location.
It's not strictly needed and could be a different directory.
