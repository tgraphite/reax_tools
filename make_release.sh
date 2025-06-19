if [[ ! -d "reax_tools_linux" ]]; then
    mkdir reax_tools_linux
    cp -r lib reax_tools_linux/
    cp doc/tutorial_1_energetic_materials.md reax_tools_linux/
    cp doc/tutorial_1_energetic_materials.pdf reax_tools_linux/
fi

if [[ ! -f "build/reax_tools" ]]; then
    rm -rf build
    cmake -B build
    cmake --build build
fi

yes | cp build/reax_tools reax_tools_linux/
tar -cjvf reax_tools_linux.tar.gz reax_tools_linux

# rm -rf reax_tools_linux