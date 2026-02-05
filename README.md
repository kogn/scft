# SCFT Solver / 自洽场理论求解器

[English](#english) | [中文](#chinese)

---

<a name="english"></a>
## English

### Overview

This is a high-performance **Self-Consistent Field Theory (SCFT)** solver designed for computational polymer physics simulations. The solver can handle complex polymer systems with spatial and orientational degrees of freedom, supporting both 1D and 2D spatial domains.

### Key Features

- **Multi-dimensional Support**: Configurable 1D and 2D spatial domains
- **Advanced Solvers**: Implementation of multiple iteration methods:
  - Picard iteration
  - Steepest Descent
  - Anderson mixing
- **Fourier Space Operations**: Leverages FFTW3 for efficient spatial transformations
- **SO3 Transformations**: Support for orientation-dependent calculations using spherical harmonics
- **High-Performance Computing**: 
  - OpenMP parallelization
  - MKL or OpenBLAS support for optimized linear algebra
  - Multi-threaded FFTW operations
- **Mathematica Integration**: WSTP (Wolfram Symbolic Transfer Protocol) support for advanced analysis

### Dependencies

Required libraries:
- **FFTW3**: Fast Fourier Transform library
- **LAPACKE**: Linear algebra routines
- **OpenBLAS** or **Intel MKL**: High-performance BLAS/LAPACK implementations
- **Mathematica** (optional): For WSTP integration

System requirements:
- C/C++ compiler with C++11 support (gcc/g++)
- OpenMP support for parallelization

### Building the Project

1. **Configure build options** by editing the `Makefile`:
   ```makefile
   DEBUG = 0        # Set to 1 for debug build
   MKL = 1          # Set to 1 for Intel MKL, 0 for standard FFTW3
   DIM = 2          # Set dimension: 1 or 2
   ```

2. **Set library paths** (update version numbers to match your installation):
   ```makefile
   INTEL_DIR = /usr/local/intel/
   MATHEMATICA_DIR = /usr/local/Wolfram/Mathematica/<version>/
   ```

3. **Build all executables**:
   ```bash
   make all
   ```

   This creates:
   - `main`: Main SCFT solver
   - `test_solver`: Solver unit tests
   - `test_trans`: Transform tests
   - `test_config`: Configuration reader tests

4. **Clean build**:
   ```bash
   make clean
   ```

### Configuration File Format

The solver uses text-based configuration files (examples in `cfg/` directory):

```
alpha = 2                  # Polymer parameter α
beta = 2                   # Polymer parameter β
kappa = 10                 # Bending rigidity
tau = 5                    # Nematic coupling parameter
chiN = 6                   # Flory-Huggins χN parameter
domain0 = 2                # Domain size in dimension 0
domain1 = 2                # Domain size in dimension 1
Steps_on_chain = 100       # Discretization steps along polymer chain
Band_width = 8             # Spectral band width
Grid_Size_x = 32           # Grid points in x-direction
Grid_Size_y = 32           # Grid points in y-direction
Output_dir = data/         # Output directory
Input_dir = data/          # Input directory
Output_filename = result   # Output file name
Input_filename = init      # Input file name
Param_filename = params    # Parameter output file
Max_steps = 100            # Maximum iteration steps
```

### Usage

**Basic usage**:
```bash
./main cfg/2D.conf
```

The solver will:
1. Read configuration from the specified file
2. Initialize the field from input file (if specified)
3. Iterate to find self-consistent solution
4. Output results to the specified data directory

**Configuration examples**:
- `cfg/1D.conf`: 1D simulation configuration
- `cfg/2D.conf`: 2D simulation configuration

### Output Files

The solver generates several output files:
- Field data files: Spatial distribution of fields
- Parameter files: System parameters and convergence information
- PDF files: Probability distribution functions

### Project Structure

```
scft/
├── src/           # Source code
│   ├── main.cpp       # Main entry point
│   ├── solver.cpp/h   # SCFT solver implementation
│   ├── iterator.cpp/hpp  # Iteration methods
│   ├── transform.cpp/h   # Fourier and SO3 transforms
│   ├── Config.cpp/h   # Configuration file parser
│   └── matrix.c/h     # Matrix operations
├── cfg/           # Configuration files
├── Makefile       # Build configuration
└── README.md      # This file
```

### Algorithm

The solver implements the standard SCFT algorithm:
1. Initialize field configuration
2. Solve modified diffusion equation for chain propagators
3. Calculate density distributions
4. Update fields based on incompressibility and interaction constraints
5. Iterate until convergence (self-consistency)

### References

For theoretical background on Self-Consistent Field Theory:
- Fredrickson, G. H. "The Equilibrium Theory of Inhomogeneous Polymers" (2006)
- Matsen, M. W. "The standard Gaussian model for block copolymer melts" (2002)

---

<a name="chinese"></a>
## 中文

### 项目概述

这是一个高性能的**自洽场理论（SCFT）求解器**，专为聚合物物理学的计算模拟设计。该求解器可以处理具有空间和取向自由度的复杂聚合物体系，支持一维和二维空间域。

### 主要特性

- **多维支持**：可配置的一维和二维空间域
- **高级求解器**：实现多种迭代方法：
  - Picard迭代
  - 最速下降法
  - Anderson混合方法
- **傅里叶空间操作**：利用FFTW3进行高效的空间变换
- **SO3变换**：使用球谐函数支持取向相关计算
- **高性能计算**：
  - OpenMP并行化
  - 支持MKL或OpenBLAS以优化线性代数运算
  - 多线程FFTW操作
- **Mathematica集成**：支持WSTP（Wolfram符号传输协议）以进行高级分析

### 依赖项

必需的库：
- **FFTW3**：快速傅里叶变换库
- **LAPACKE**：线性代数例程
- **OpenBLAS** 或 **Intel MKL**：高性能BLAS/LAPACK实现
- **Mathematica**（可选）：用于WSTP集成

系统要求：
- 支持C++11的C/C++编译器（gcc/g++）
- 支持OpenMP的并行化

### 构建项目

1. **配置构建选项**，编辑`Makefile`：
   ```makefile
   DEBUG = 0        # 设置为1进行调试构建
   MKL = 1          # 设置为1使用Intel MKL，0使用标准FFTW3
   DIM = 2          # 设置维度：1或2
   ```

2. **设置库路径**（更新版本号以匹配您的安装）：
   ```makefile
   INTEL_DIR = /usr/local/intel/
   MATHEMATICA_DIR = /usr/local/Wolfram/Mathematica/<version>/
   ```

3. **构建所有可执行文件**：
   ```bash
   make all
   ```

   这将创建：
   - `main`：主SCFT求解器
   - `test_solver`：求解器单元测试
   - `test_trans`：变换测试
   - `test_config`：配置读取器测试

4. **清理构建**：
   ```bash
   make clean
   ```

### 配置文件格式

求解器使用基于文本的配置文件（示例在`cfg/`目录中）：

```
alpha = 2                  # 聚合物参数α
beta = 2                   # 聚合物参数β
kappa = 10                 # 弯曲刚度
tau = 5                    # 向列耦合参数
chiN = 6                   # Flory-Huggins χN参数
domain0 = 2                # 维度0的域大小
domain1 = 2                # 维度1的域大小
Steps_on_chain = 100       # 沿聚合物链的离散化步数
Band_width = 8             # 谱带宽
Grid_Size_x = 32           # x方向的网格点数
Grid_Size_y = 32           # y方向的网格点数
Output_dir = data/         # 输出目录
Input_dir = data/          # 输入目录
Output_filename = result   # 输出文件名
Input_filename = init      # 输入文件名
Param_filename = params    # 参数输出文件
Max_steps = 100            # 最大迭代步数
```

### 使用方法

**基本用法**：
```bash
./main cfg/2D.conf
```

求解器将：
1. 从指定文件读取配置
2. 从输入文件初始化场（如果指定）
3. 迭代以找到自洽解
4. 将结果输出到指定的数据目录

**配置示例**：
- `cfg/1D.conf`：一维模拟配置
- `cfg/2D.conf`：二维模拟配置

### 输出文件

求解器生成几个输出文件：
- 场数据文件：场的空间分布
- 参数文件：系统参数和收敛信息
- PDF文件：概率分布函数

### 项目结构

```
scft/
├── src/           # 源代码
│   ├── main.cpp       # 主入口点
│   ├── solver.cpp/h   # SCFT求解器实现
│   ├── iterator.cpp/hpp  # 迭代方法
│   ├── transform.cpp/h   # 傅里叶和SO3变换
│   ├── Config.cpp/h   # 配置文件解析器
│   └── matrix.c/h     # 矩阵操作
├── cfg/           # 配置文件
├── Makefile       # 构建配置
└── README.md      # 本文件
```

### 算法

求解器实现标准的SCFT算法：
1. 初始化场配置
2. 求解链传播子的修正扩散方程
3. 计算密度分布
4. 基于不可压缩性和相互作用约束更新场
5. 迭代直到收敛（自洽）

### 参考文献

关于自洽场理论的理论背景：
- Fredrickson, G. H. "The Equilibrium Theory of Inhomogeneous Polymers" (2006)
- Matsen, M. W. "The standard Gaussian model for block copolymer melts" (2002)
