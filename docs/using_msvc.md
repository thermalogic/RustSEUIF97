# Using MSVC to build the lib

##  Building 

```bash
cargo build -r --features cdecl
```

- stdcall: Windows API functions(MSVC 64bit)

```bash
cargo build -r --features stdcall
```

- stdcall: Windows API functions(MSVC 32bit)

```bash
cargo build -r  --target=i686-pc-windows-msvc --features stdcall
```

## Building test

* cmake

cd build
```bash
cmake ..
cmake --build . --config release
``
