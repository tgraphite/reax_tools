if [[ -d "reax_tools_linux" ]]; then
    rm -rf reax_tools_linux
fi

if [[ -d "reax_tools_windows" ]]; then
    rm -rf reax_tools_windows
fi

mkdir reax_tools_linux

cp build/reax_tools reax_tools_linux/
cp -r lib reax_tools_linux/
cp doc/tutorial_1_energetic_materials.md reax_tools_linux/
cp doc/tutorial_1_energetic_materials.pdf reax_tools_linux/
cp plot/* reax_tools_linux/

tar -cjvf reax_tools_linux.tar.gz reax_tools_linux

rm -rf reax_tools_linux