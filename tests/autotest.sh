time reax_tools -f cp2k.xyz -me C -r 1.25
time reax_tools -f FeHCONSF.xyz -t X,H,C,O,N,S,F -r 1.2 -nt 4 -me C -rc
time reax_tools -f CFOH.lammpstrj -t C,F,O,H -r 1.2
time reax_tools -f CNOH.lammpstrj -t C,N,O,H -r 1.2
time reax_tools -f CHO.lammpstrj -t C,H,O -r 1.2 -me C

for dir in cp2k_reax_tools FeHCONSF_reax_tools CFOH_reax_tools CNOH_reax_tools CHO_reax_tools; do
    cd $dir
    dot -Tpng reaction_flow.dot -o reaction_flow.png
    cd ..
done
