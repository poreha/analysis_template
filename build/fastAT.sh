#!/bin/bash
# Builds AT template prog, performs analysis and opens rootbrowse
cmake ..
make -j2
./analyse ../lists/AgAg_1_58A_GeV.list
rootbrowse output.root
