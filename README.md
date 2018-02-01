# simpletools
Simpletools: Handy command line tools for ntuple manipulation and analysis  \
&copy; Conor Fitzpatrick 2007-2015

If you find these tools useful either in whole or in part, please cite the documentation:

>simpletools: Handy command line tools for ntuple manipulation and analysis  \
>[LHCb Internal Note 2009-029](https://cds.cern.ch/record/1223091)  \
>Conor Fitzpatrick, University of Edinburgh

The most recent version of simpletools, including the above documentation will always be available from  \
`/afs/cern.ch/user/c/cofitzpa/public/simpletools/`

Bug reports, feature requests, comments and patches are always welcome.  \
conor.fitzpatrick@cern.ch

## Installation

Building and installing are done with CMake.
If you want to install simpletools system-wide, run:

    cmake .
    make
    sudo make install

If you don't have administrator privileges and/or want to install locally, use `-DCMAKE_INSTALL_PREFIX` to change the destination, *e.g.*:

    cmake . -DCMAKE_INSTALL_PREFIX=~/.local/
    make
    make install

### Uninstallation

After running `make install`, a file called `install_manifest.txt` will be created, which lists the destinations of all installed binaries and libraries.
To uninstall, you can run:

    make uninstall

which will delete every file listed in `install_manifest.txt`.

