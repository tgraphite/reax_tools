#!/bin/bash

cd ../build
rm -rf *
cmake ..
cmake --build . -j8

cd ../examples

../build/reax_tools -f cp2k.xyz -me O
../build/reax_tools -f small.lammpstrj -t Fe,O,H -me O
../build/reax_tools -f drill.lammpstrj -t X,H,C,O,N,S,F -me C --dump
../build/reax_tools -f polymer.lammpstrj -t C,F,O,H -me C -rc

rm -rf *.png

dot -Tpng cp2k.dot > cp2k.png
dot -Tpng small.dot > small.png
dot -Tpng drill.dot > drill.png
dot -Tpng polymer.dot > polymer.png

