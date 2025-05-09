cd .. 
# rm -rf build
# cmake -B build 
cmake --build build -j8
yes | cp build/reax_tools reax_tools
cd tests

time ../reax_tools/reax_tools -f cp2k.xyz -me C -r 1.25
time ../reax_tools/reax_tools -f FeHCONSF.xyz -t Fe,H,C,O,N,S,F -r 1.2 -nt 4 -me C -rc
time ../reax_tools/reax_tools -f CFOH.lammpstrj -t C,F,O,H -r 1.2 --dump

for dir in cp2k_reax_tools FeHCONSF_reax_tools CFOH_reax_tools; do
    cd $dir
    dot -Tpng reaction_flow.dot -o reaction_flow.png
    cd ..
done
