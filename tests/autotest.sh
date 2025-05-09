cd .. 
# rm -rf build
# cmake -B build 
cmake --build build -j8
yes | cp build/reax_tools reax_tools
cd tests

time ../reax_tools/reax_tools -f cp2k.xyz -me C
time ../reax_tools/reax_tools -f drill.lammpstrj -t Fe,H,C,O,N,S,F -r 1.2 -nt 4 -me C -rc
time ../reax_tools/reax_tools -f polymer.lammpstrj -t C,F,O,H -r 1.15
time ../reax_tools/reax_tools -f small.lammpstrj -t Fe,O,H -r 1.15

