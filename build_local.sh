# Remove old build files
rm -f test/reax_tools

SRC_FILES=$(find cpp -name '*.cpp')
SRC_FILES="$SRC_FILES cpp/fmt/format.cc cpp/fmt/os.cc cpp/universe_parallel.cc cpp/local_main.cc cpp/rdkit_utils.cc"
SRC_FILES=${SRC_FILES/cpp\/wasm_main.cpp/}

echo $SRC_FILES
echo "=================="

g++ $SRC_FILES -o bin/reax_tools \
  -Icpp -Icpp/fmt \
  -O2 \
  -std=c++17 \
  -I/opt/cpp/rdkit/include/rdkit -L/opt/cpp/rdkit/lib \
  -I/opt/cpp/boost_1_88_0/include/boost -L/opt/cpp/boost_1_88_0/lib \
  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitSmilesParse -lRDKitFileParsers -lRDKitRDGeometryLib -lRDKitRDStreams -lRDKitDistGeomHelpers -lRDKitMolDraw2D 

