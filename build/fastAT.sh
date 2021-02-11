#!/bin/bash
# Builds AT template prog, performs analysis and opens rootbrowse
make -j2
./analyse ../lists/AuAu_1_23AGeV_GEN9_v4.list
rootbrowse output.root
