# rm -rf ../../build_linux/*
cd ../../build_linux
# cmake ..
cmake --build . -j 8
cd ../examples/linux
../../build_linux/reax_tools -f cp2k.xyz -me C
# ../../build_linux/reax_tools -f drill.lammpstrj -t Fe,H,C,O,N,S,F -r 1.2 -nt 4 -me C -rc --dump
# ../../build_linux/reax_tools -f polymer.lammpstrj -t C,F,O,H -r 1.15
# ../../build_linux/reax_tools -f small.lammpstrj -t Fe,O,H -r 1.15