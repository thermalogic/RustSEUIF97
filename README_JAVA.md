# Java Using the Shared Library building from Rust

This document describes two methods for Java to call the Rust shared library, with support for **Windows, Linux, macOS** platforms.

1. **JNA (Java Native Access)** - Simple and widely compatible
2. **Panama FFI** - Built-in from Java 22+, more modern approach

## Prerequisites

### 1. Build Rust Library

| Platform | Command | Output Library |
|----------|---------|----------------|
| **Windows** | `cargo build -r --features cdecl --target x86_64-pc-windows-msvc` | `target/x86_64-pc-windows-msvc/release/seuif97.dll` |
| **Linux** | `cargo build -r --features cdecl --target x86_64-unknown-linux-gnu` | `target/x86_64-unknown-linux-gnu/release/libseuif97.so` |
| **macOS** | `cargo build -r --features cdecl --target x86_64-apple-darwin` | `target/x86_64-apple-darwin/release/libseuif97.dylib` |

**Shorter command** (uses default target for your platform):
```bash
cargo build -r --features cdecl
```

### 2. Java Environment Requirements

| Method | Minimum Java Version | Additional Dependencies |
|--------|---------------------|------------------------|
| JNA | Java 5+ | jna.jar |
| Panama FFI | Java 22+ | None (built-in) |

---

## Method 1: JNA (Recommended)

JNA is the most common way to call native libraries from Java. It's simple and has good cross-platform compatibility.

### Steps

1. **Download JNA**

   Download `jna.jar` from [Maven Central](https://mvnrepository.com/artifact/net.java.dev.jna/jna) and place it in the `demo_using_lib` directory.

2. **Compile and Run**

**Windows:**
```bash
cd demo_using_lib

# Compile
javac -cp jna.jar demo_jni_cdecl.java

# Run
java -cp jna.jar;. -Djava.library.path=../target/x86_64-pc-windows-msvc/release demo_jni_cdecl
```

**Linux:**
```bash
cd demo_using_lib

# Compile
javac -cp jna.jar demo_jni_cdecl.java

# Run
java -cp jna.jar:. -Djava.library.path=../target/x86_64-unknown-linux-gnu/release demo_jni_cdecl
```

### Code Example

```java
import com.sun.jna.Library;
import com.sun.jna.Native;

public interface SEUIF97 extends Library {
    SEUIF97 INSTANCE = Native.load("seuif97", SEUIF97.class);
    
    double pt(double p, double t, int o_id);
 }

public class demo_jni_cdecl {
    public static void main(String[] args) {
        double h = SEUIF97.INSTANCE.pt(16.0, 540.0, 4);  // enthalpy
        System.out.printf("Enthalpy: %.3f kJ/kg%n", h);
    }
}
```

---

## Method 2: Panama FFI (Java 22+)

Panama FFI is the Foreign Function & Memory API built into Java 22+. No additional dependencies required.

### Steps

**Windows:**
```bash
cd demo_using_lib

# Compile
javac seuif97_panama_simple.java

# Run
java -Djava.library.path=../target/x86_64-pc-windows-msvc/release seuif97_panama_simple
```

**Linux:**
```bash
cd demo_using_lib

# Compile
javac seuif97_panama_simple.java

# Run
java -Djava.library.path=../target/x86_64-unknown-linux-gnu/release seuif97_panama_simple
```

### Code Example

```java
import java.lang.foreign.*;
import java.lang.invoke.MethodHandle;

public class seuif97_panama_simple {
    static {
        System.loadLibrary("seuif97");
    }
    
    private static final Linker LINKER = Linker.nativeLinker();
    private static final SymbolLookup LOOKUP = SymbolLookup.loaderLookup();
    
    private static final MethodHandle PT;
    
    static {
        try {
            MemorySegment ptSeg = LOOKUP.find("pt").orElseThrow();
            FunctionDescriptor desc = FunctionDescriptor.of(
                ValueLayout.JAVA_DOUBLE,
                ValueLayout.JAVA_DOUBLE,
                ValueLayout.JAVA_DOUBLE,
                ValueLayout.JAVA_INT
            );
            PT = LINKER.downcallHandle(ptSeg, desc);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
    
    public static double pt(double p, double t, int o_id) throws Throwable {
        return (double) PT.invokeExact(p, t, o_id);
    }
}
```

## Comparison Table

| Aspect | JNA | Panama FFI |
|--------|-----|------------|
| **Java Version** | 5+ | 22+ |
| **Dependencies** | jna.jar required | None (built-in) |
| **Code Complexity** | Low (declare interface only) | Slightly higher (manual configuration) |
| **Performance** | Slightly slower (JNA wrapper) | Better (direct FFI) |
| **Learning Curve** | Gentle | Steeper |
| **Compatibility** | Excellent (all platforms) | Limited (Java 22+ only) |
| **Maven Publishing** | Simple | Simple (no dependencies) |

## Recommendations

| Scenario | Recommended Method |
|----------|-------------------|
| Production environment, wide compatibility needed | **JNA** |
| New project, using Java 22+ | **Panama FFI** |
| Need to publish to Maven | **Panama FFI** (no dependencies) |
| Rapid prototyping | **JNA** (simpler) |

## Run Scripts

### Windows Scripts

**JNA (`run_jni.bat`)**
```batch
@echo off
set LIB_PATH=../target/x86_64-pc-windows-msvc/release
javac -cp jna.jar demo_jni_cdecl.java
java -cp jna.jar;. -Djava.library.path="%LIB_PATH%" demo_jni_cdecl
pause
```

**Panama FFI (`run_panama_simple.bat`)**
```batch
@echo off
set LIB_PATH=../target/x86_64-pc-windows-msvc/release
javac seuif97_panama_simple.java
java -Djava.library.path="%LIB_PATH%" seuif97_panama_simple
pause
```

### Linux Scripts

**JNA (`run_jni.sh`)**
```bash
#!/bin/bash
LIB_PATH="../target/x86_64-unknown-linux-gnu/release"
javac -cp jna.jar demo_jni_cdecl.java
java -cp jna.jar:. -Djava.library.path="$LIB_PATH" demo_jni_cdecl
```

**Panama FFI (`run_panama_simple.sh`)**
```bash
#!/bin/bash
LIB_PATH="../target/x86_64-unknown-linux-gnu/release"
javac seuif97_panama_simple.java
java -Djava.library.path="$LIB_PATH" seuif97_panama_simple
```
