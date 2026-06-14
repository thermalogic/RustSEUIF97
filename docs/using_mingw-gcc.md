# Windows 下使用 MinGW-GCC 编译 Rust C 绑定共享库指南

## 1. 环境准备

### 1.1 安装 Rust 和 GNU Target

```bash
# 安装 GNU target
rustup target add x86_64-pc-windows-gnu
```

### 1.2 安装 MinGW-w64

确保 MinGW-w64 已安装，并将 `bin` 目录添加到系统 PATH 环境变量中。

验证安装：

```bash
gcc --version
```

## 2. 配置 Cargo

在项目根目录创建 `.cargo/config.toml`：

```toml
[target.x86_64-pc-windows-gnu]
linker = "x86_64-w64-mingw32-gcc"
ar = "x86_64-w64-mingw32-gcc-ar"
```

## 3. 编译 Rust 共享库

### 3.1 构建命令

```bash
cargo build -r --features cdecl --target x86_64-pc-windows-gnu
```

### 3.2 输出文件

编译完成后，输出文件位于 `target/x86_64-pc-windows-gnu/release/`：

- `seuif97.dll` - 动态链接库
- `libseuif97.dll.a` - 导入库（用于链接）

**注意**：Rust 的 `cdylib` 类型在 Windows GNU target 上生成的导入库文件名为 `libseuif97.dll.a`。

## 4. 使用 MinGW-GCC 链接共享库

### 4.1 使用导入库（推荐）

Rust 编译后会生成 `libseuif97.dll.a` 导入库，直接使用 `-l` 参数链接：

```bash
gcc -o myapp.exe myapp.c -L target/x86_64-pc-windows-gnu/release -lseuif97
```

**说明**：`-lseuif97` 会自动查找 `libseuif97.dll.a`。

## 5. 使用示例

### 5.1 C 代码示例 (`test_seuif97.c`)

```c
#include <stdio.h>

// 声明 Rust 库中的函数（根据实际导出的函数调整）
extern double h_pt(double p, double t);

int main() {
    double p = 16.0;  // MPa
    double t = 500.0; // °C
    
    double h = h_pt(p, t);
    printf("Enthalpy at p=%.1f MPa, t=%.1f °C: h = %.2f kJ/kg\n", p, t, h);
    
    return 0;
}
```

### 5.2 编译命令

```bash
gcc -o test_seuif97.exe test_seuif97.c -L target/x86_64-pc-windows-gnu/release -lseuif97
```

### 5.3 运行

```bash
# 确保 seuif97.dll 在可执行文件同一目录或 PATH 中
copy target\x86_64-pc-windows-gnu\release\seuif97.dll .

.\test_seuif97.exe
```

## 6. Makefile 示例

```makefile
CC = gcc
CFLAGS = -Wall -O2
LDFLAGS = -L target/x86_64-pc-windows-gnu/release -lseuif97
TARGET = test_seuif97.exe
SRC = test_seuif97.c

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

run: $(TARGET)
	.\$(TARGET)

clean:
	del /Q $(TARGET) *.dll *.a *.def 2>nul
```

