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
vertices, faces, tetrahedra, is_apex = cgal_alpha.compute_alpha_complex_from_atoms(
    positions, radii, alpha_value, probe_radius, filter_type
)

print(f"Extracted {vertices.shape[0]} vertices ({is_apex.sum()} wire apexes)")
print(f"Extracted {faces.shape[0]} faces")
print(f"Extracted {tetrahedra.shape[0]} interior tetrahedra")
```

## Output

- `vertices` (N, 3) float32: one vertex per atom on a face or tetrahedron, then one apex vertex per wire triangle.
- `faces` (M, 3) int32: the facets whose classification passes `filter`, then the wire triangles.
- `tetrahedra` (K, 4) int32: the INTERIOR cells, indexing atom vertices.
- `is_apex` (N,) bool: true for the apex vertices, which are always the last ones.

A linking wire is a SINGULAR edge whose two ends lie on faces of two different connected components of the face mesh. Each one becomes the triangle (u, v, apex), with the apex at the neck radius rho over the edge midpoint, rho being the radius of the circle where the two probe-inflated balls (radius r + probe_radius) meet. The apex direction is `normalize(cross(v - u, helper))`, with `helper` the x axis unless `|(v - u).x| >= 0.9 |v - u|`, in which case it is the y axis.

## License
Provided as-is. Check CGAL's licensing requirements if using in a commercial application.
