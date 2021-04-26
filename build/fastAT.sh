#!/bin/bash
# Builds AT template prog, performs analysis and opens rootbrowse
# AgAg_1_58A_GeV - is MC-generator + Geant3
cmake ..
make -j2
./analyse ../lists/AgAg_1_58A_GeV.list
rootbrowse output.root
