# SEUIF97 C Binding

SEUIF97 provides C bindings through Rust's FFI (Foreign Function Interface). On Windows, you can choose between two toolchains:

- **MSVC** (Microsoft Visual C++): Native Windows toolchain, produces standard Windows DLLs
- **MinGW-GCC** (Minimalist GNU for Windows): GNU toolchain for Windows, compatible with GCC ecosystem

## Prerequisites

### MSVC Toolchain

Install [Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022) or Visual Studio with **Desktop development with C++** workload.

Verify installation:
```bash
cl
```

### MinGW-GCC Toolchain

1. Install [MinGW-w64](https://winlibs.com/) (recommended) or [MSYS2](https://www.msys2.org/)
2. Add `bin` directory to system `PATH`
3. Install Rust target:
```bash
rustup target add x86_64-pc-windows-gnu
```

Verify installation:
```bash
gcc --version
```

> **Note:** MSVC is the default toolchain on Windows. No linker configuration is needed for MSVC builds.

## Building the DLL

### Option 1: MSVC (Default, Recommended for Windows)

MSVC is the default toolchain on Windows. No additional linker configuration is required.

The library supports two calling conventions on Windows:

**cdecl** (default, compatible with most C compilers):
```bash
# Uses MSVC linker automatically (default target: x86_64-pc-windows-msvc)
cargo build -r --features cdecl
```

**stdcall** (Windows API convention, MSVC 64-bit):
```bash
cargo build -r --features stdcall
```

**stdcall for 32-bit Windows:**
```bash
cargo build -r --target=i686-pc-windows-msvc --features stdcall
```

Output location: `target/release/seuif97.dll`

### Option 2: MinGW-GCC (Requires explicit target)

When using MinGW-GCC, you must explicitly specify the target:

```bash
cargo build -r --features cdecl --target x86_64-pc-windows-gnu
```

Output location: `target/x86_64-pc-windows-gnu/release/`

Files generated:
- `seuif97.dll` - Dynamic link library
- `libseuif97.dll.a` - Import library for linking

## Using the DLL in C

### Example Code

```c
#include <stdio.h>

// Declare the function from the DLL
extern double pt2h(double p, double t);

int main() {
    double p = 16.0;   // MPa
    double t = 500.0;  // °C

    double h = pt2h(p, t);
    printf("Enthalpy at p=%.1f MPa, t=%.1f °C: h = %.2f kJ/kg\n", p, t, h);

    return 0;
}
```

### Compiling and Linking

#### With MSVC (cl)

```bash
cl demo.c /link seuif97.lib
```

Or with the DLL directly:
```bash
cl demo.c /link target\release\seuif97.dll.lib
```

#### With MinGW-GCC

```bash
gcc demo.c -L./target/x86_64-pc-windows-gnu/release -lseuif97 -o demo.exe
```

> **Note:** Ensure `seuif97.dll` is in the same directory as `demo.exe` or in a directory listed in your `PATH`.

### Using the Provided Makefile

The project includes a cross-platform `makefile`:

```makefile
ifeq ($(OS),Windows_NT)
    LIBFLAGS=-L../../target/release -lseuif97
    EXEDIR=../../target/release/
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        LIBFLAGS=-L../../target/release -Wl,-rpath=../../target/release -lseuif97 -lm
        EXEDIR=./
    endif
endif

all:
    gcc demo.c -o$(EXEDIR)demo $(LIBFLAGS)
    $(EXEDIR)demo
```

Run with:
```bash
cd demo_using_lib/demo_c
make
```

## Calling Convention Reference

| Convention | Feature Flag | Target | Use Case |
|------------|--------------|--------|----------|
| cdecl | `--features cdecl` | MSVC, MinGW-GCC | Default, compatible with most C code |
| stdcall | `--features stdcall` | MSVC 32/64-bit | Windows API compatibility |

