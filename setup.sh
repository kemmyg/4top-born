#!/bin/sh


echo 'Add binaries to the path...'
export lhapdf_config=/Users/kemmler/Programms/LHAPDF-6.4.0/BUILD/bin
export PATH=$PATH:/Users/kemmler/Programms/4top-born
echo 'Define LHALIB...'
export PATH=$PATH:/Users/kemmler/Programms/LHAPDF-6.4.0/BUILD/bin
export LD_LIBRARY_PATH=/Users/kemmler/Programms/LHAPDF-6.4.0/BUILD/lib
export PYTHONPATH=/Users/kemmler/Programms/LHAPDF-6.4.0/BUILD/lib/python2.7/site-packages/
echo 'Done!'
