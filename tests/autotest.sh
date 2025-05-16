reax_tools -f cp2k.xyz -me C -r 1.25
reax_tools -f FeHCONSF.xyz -t X,H,C,O,N,S,F -r 1.2 -nt 4 -me C -rc
reax_tools -f CFOH.lammpstrj -t C,F,O,H -r 1.2
reax_tools -f CNOH.lammpstrj -t C,N,O,H -r 1.2
reax_tools -f CHO.lammpstrj -t C,H,O -r 1.2 -me C

for dir in cp2k_rtresults FeHCONSF_rtresults CFOH_rtresults CNOH_rtresults CHO_rtresults; do
    cd $dir
    dot -Tpng reaction_flow.dot -o reaction_flow.png
    cd ..
done
