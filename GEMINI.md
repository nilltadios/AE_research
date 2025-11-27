# Project Context: Aeroelasticity & ERA Analysis

## Overview
This project contains MATLAB/Octave scripts for simulating and analyzing the aeroelastic behavior of an airfoil (pitch and plunge degrees of freedom). It focuses on identifying the onset of flutter (instability) and Limit Cycle Oscillations (LCO) using non-linear stiffness models.

A key feature is the integration of the **Eigensystem Realization Algorithm (ERA)** to reconstruct system states/outputs from time-domain data.

**Recent Update (Nov 2025):** The project has been ported to support **GPU acceleration** on Linux systems using AMD GPUs via the **OpenCL** backend for GNU Octave.

## Directory Structure

*   **`gpu_optimized/`**: The main GPU-accelerated codebases.
    *   **`octave_version/`**: **Primary Version.** Optimized for GNU Octave using the `ocl` package (OpenCL). Includes `mu_main_optimized.m`, `aeroelastic.m`, etc.
    *   **`matlab_version/`**: MATLAB-compatible variant using Parallel Computing Toolbox (`gpuArray`).
*   **`cpu_legacy/`**: Original CPU-based analysis scripts (`mu_main.m`, etc.). Kept for reference, validation, and machines without GPU support.
*   **`ERA_Embedded/`**: Independent implementation targeting embedded systems (self-contained).
*   **`tests/`**: Validation suites.
    *   **`ERA/`**: Tests for the ERA algorithm.

## GPU Acceleration (`gpu_optimized/octave_version/`)

### Key Features
*   **OpenCL Backend:** Utilizes the `ocl` Octave package to interface with AMD/Intel GPUs via Mesa/Rusticl drivers.
*   **Performance Optimization:** 
    *   Vectorized critical speed calculations.
    *   Reduced redundant object creation in ODE solvers.
    *   Optimized Hankel matrix construction.
    *   **Randomized SVD (`rsvd_gpu`)**: Support for randomized Singular Value Decomposition.
*   **Robust Hybrid Pipeline:** Automatically detects available compute backends.

### Usage

1.  **Prerequisites:**
    *   GNU Octave
    *   `ocl` package (`pkg load ocl`)
    *   `signal` package (`pkg load signal`)
    *   OpenCL drivers

2.  **Execution:**
    Run the optimized script directly from the terminal:
    ```bash
    octave --eval "pkg load signal ocl; cd gpu_optimized/octave_version; mu_main_optimized"
    ```

3.  **Configuration:**
    *   In `mu_main_optimized.m`, set `PERFORMANCE_MODE = true` for faster execution.
    *   Adjust `era_params` for trade-off between speed and accuracy.

## Standard Components (CPU)

To run the legacy CPU version, navigate to `cpu_legacy/` and run `mu_main.m`.

### Main Files
*   **`mu_main.m`**: The legacy CPU entry point.
*   **`aeroelastic.m`**: System ODEs.

## Notes
*   **Stability:** The GPU pipeline currently relies on CPU fallback for specific SVD steps due to driver limitations with OpenCL `float matrix` types on consumer AMD cards.
*   **Compatibility:** The `findpeaks` implementation has been patched to handle negative values correctly in Octave.