# CGAL Alpha Complex Python Bindings

A fast and efficient Python wrapper for computing 3D alpha complexes (volumes and boundary surfaces) from atom positions and radii. This library utilizes the highly optimized `CGAL::Fixed_alpha_shape_3` and is exposed to Python via pybind11. Thin necks of the surface, carried by edges that no facet contains, are added as wire triangles so that the components they join stay connected.

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

`compute_alpha_complex_from_atoms` returns `(vertices, faces, tetrahedra, is_apex)`:

- `vertices` (N, 3) float32: the atom centres used by a face, a tetrahedron or a wire triangle, followed by one apex vertex per wire triangle.
- `faces` (M, 3) int32: the facets whose classification passes `filter`, followed by the wire triangles.
- `tetrahedra` (K, 4) int32: the INTERIOR cells, indexing atom vertices only.
- `is_apex` (N,) bool: marks the apex vertices, the last `is_apex.sum()` rows of `vertices`.

The surface is not tufted: each atom has a single vertex, shared by faces, wire triangles and tetrahedra.

## Wire triangles

A wire is a SINGULAR edge of the alpha complex, an edge that no facet of the complex contains. It carries a neck of the solvent-accessible surface: the circle of radius rho where the two balls of radius r + probe_radius meet. Wires that link two connected components of the face mesh are kept:

- directly, when both ends lie on faces of different components;
- through a chain, when the wire touches a lone atom (on no face) and its cluster of lone atoms, joined by wires, reaches two or more components. Wires ending at a lone atom that has no other wire are removed first, repeatedly.

Each kept wire (u, v) becomes the triangle (u, v, apex), with

```
apex = (u + v) / 2 + rho * normalize(cross(v - u, helper))
```

where `helper` is the x axis, or the y axis when `|(v - u).x| >= 0.9 |v - u|`.

## License
Provided as-is. Check CGAL's licensing requirements if using in a commercial application.
