rm -f web/wasm_main.*

# Collect all source files (excluding fmt.cc, only compile format.cc and os.cc)
SRC_FILES=$(find cpp -name '*.cpp')
SRC_FILES="$SRC_FILES cpp/fmt/format.cc cpp/fmt/os.cc"
# SRC_FILES=${SRC_FILES/}

# 编译
emcc $SRC_FILES -o web/wasm_main.js \
  -O3 \
  -Icpp -Icpp/fmt \
  -s WASM=1 \
  -s "EXPORTED_RUNTIME_METHODS=['ccall','FS','lengthBytesUTF8','stringToUTF8','setValue']" \
  -s "EXPORTED_FUNCTIONS=['_cpp_main','_malloc','_free']" \
  -s FORCE_FILESYSTEM=1 \
  -s INITIAL_MEMORY=512MB \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s ASSERTIONS=1 \
  -s DISABLE_EXCEPTION_CATCHING=0 \
  -s INVOKE_RUN=0 \
  -D WASM_MODE=1 \