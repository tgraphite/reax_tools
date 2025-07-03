SRC_FILES=$(find cpp -name '*.cpp')
SRC_FILES="$SRC_FILES cpp/fmt/format.cc cpp/fmt/os.cc cpp/universe_parallel.cc cpp/local_main.cc cpp/rdkit_utils.cc"
SRC_FILES=${SRC_FILES/cpp\/wasm_main.cpp/}

echo $SRC_FILES
echo "=================="

g++ $SRC_FILES -o bin/reax_tools \
  -Icpp -Icpp/fmt \
  -O3 \
  -std=c++17 \
  -DENABLE_RDKIT \
  -Wl,-rpath,'$ORIGIN/lib' \
  -I/opt/cpp/rdkit/include/rdkit -L/opt/cpp/rdkit/lib \
  -I/opt/cpp/boost_1_88_0/include/boost -L/opt/cpp/boost_1_88_0/lib \
  -lRDKitGraphMol -lRDKitRDGeneral -lRDKitSmilesParse -lRDKitFileParsers -lRDKitRDGeometryLib -lRDKitRDStreams -lRDKitDistGeomHelpers -lRDKitMolDraw2D 

# Find and copy all non-system .so dependencies to bin/
DEPS=$(ldd bin/reax_tools | grep -E '/opt/cpp|/home|/usr/local' | awk '{print $3}')
for dep in $DEPS; do
  cp -u "$dep" bin/lib
done


