# CGAL Alpha Complex Python Bindings

A fast and efficient Python wrapper for computing 3D alpha complexes (volumes and boundary surfaces) from atom positions and radii. This library utilizes the highly optimized `CGAL::Fixed_alpha_shape_3` and is exposed to Python via pybind11.

## Dependencies

- **CMake** (>= 3.16)
- **C++20** compatible compiler
- **Python 3** and **pybind11**
- **CGAL** (including GMP and MPFR)

## Installation

You can build this module using CMake:

```bash
mkdir build && cd build
cmake ..
make
```

Ensure the compiled `.so` (or `.pyd` on Windows / `.dylib` on macOS) is located in your python path.

## Usage

```python
import numpy as np
import cgal_alpha

# Setup mock data (N atoms)
N = 100
positions = np.random.rand(N, 3).astype(np.float32) * 10.0
radii = np.ones(N, dtype=np.float32) * 1.5

alpha_value = 0.0          # The alpha threshold (squared radius for Fixed_alpha_shape_3)
probe_radius = 1.4         # Used to expand atomic radii (e.g., solvent probe radius)
filter_type = "singular+regular" 

# Compute the alpha complex
vertices, faces, tetrahedra = cgal_alpha.compute_alpha_complex_from_atoms(
    positions, radii, alpha_value, probe_radius, filter_type
)

print(f"Extracted {vertices.shape[0]} vertices")
print(f"Extracted {faces.shape[0]} faces")
print(f"Extracted {tetrahedra.shape[0]} interior tetrahedra")
```

## License
Provided as-is. Check CGAL's licensing requirements if using in a commercial application.
