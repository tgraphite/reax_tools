cd ../../build
rm -rf *
cmake ..
cmake --build . -j 4
cd ../examples/linux