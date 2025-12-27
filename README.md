A high-performance, CPU-only Minecraft-style renderer. This project is a "just for fun" technicaldemonstration of the relevance of the software renderer in 2025, bypassing the GPU entirely to render 3D scenes using pure C++, SIMD, and multi-threading.

## 🚀 Performance Benchmarks

The renderer scales across different CPU architectures by leveraging specific SIMD instruction sets.

### 1
* **CPU:** AMD Ryzen 5 2700 (8C/16T)
* **Instruction Set:** SSE4.2 (128-bit)
* **Resolution:** 1920x1080 | 15 Chunks
* **Result:** 50–60 FPS
* *Observation:* High throughput achieved through massive multi-threading despite narrower SIMD registers.

### 2
* **CPU:** Intel Core i5-8350U (4C/8T)
* **Instruction Set:** AVX2 (256-bit)
* **Resolution:** 1920x1080 | 15 Chunks
* **Result:** 30–40 FPS
* *Observation:* The use of AVX2 significantly compensates for the lower core count, processing twice as many pixels per clock cycle compared to SSE.

## 📦 Requirements & Building
The project uses **vcpkg** for dependency management.
### Dependencies  
* C++20
* CMake
* SDL2
* xsimd
* glad
* glm
* Stb
* OpenMP

## Building
If you have all dependencies installed in your system (via vcpkg, brew, or apt), it should "just work".

---

![Preview](/photo.jpg)
![Preview](/phono1.png)
