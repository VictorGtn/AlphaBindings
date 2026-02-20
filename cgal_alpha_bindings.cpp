// CGAL Alpha Complex Python Bindings
// pybind11 wrapper using CGAL's Fixed_alpha_shape_3 for optimal performance

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>

#include <vector>
#include <array>
#include <map>
#include <stdexcept>

namespace py = pybind11;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_vertex_base_3<K> Vbb;
typedef CGAL::Fixed_alpha_shape_vertex_base_3<K, Vbb> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K> Rcb;
typedef CGAL::Fixed_alpha_shape_cell_base_3<K, Rcb> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<K, Tds> Triangulation;
typedef CGAL::Fixed_alpha_shape_3<Triangulation> Fixed_alpha_shape;

typedef Triangulation::Weighted_point Weighted_point;
typedef Triangulation::Bare_point Point_3;
typedef Fixed_alpha_shape::Vertex_handle Vertex_handle;

// Returns (vertices, faces, tetrahedra).
// faces:        boundary triangles (REGULAR facets, i.e. between INTERIOR and EXTERIOR)
// tetrahedra:   volume cells classified as INTERIOR
std::tuple<py::array_t<float>, py::array_t<int32_t>, py::array_t<int32_t>>
extract_alpha_volume(Fixed_alpha_shape& A) {
    std::map<Vertex_handle, int> vertex_index;
    std::vector<Vertex_handle> vertices;

    auto get_or_insert = [&](Vertex_handle vh) -> int {
        auto it = vertex_index.find(vh);
        if (it != vertex_index.end()) return it->second;
        int idx = vertices.size();
        vertex_index[vh] = idx;
        vertices.push_back(vh);
        return idx;
    };

    std::vector<std::array<int, 3>> faces;
    std::vector<std::array<int, 4>> tets;

    // Boundary triangles: REGULAR facets
    for (auto fit = A.finite_facets_begin(); fit != A.finite_facets_end(); ++fit) {
        if (A.classify(*fit) != Fixed_alpha_shape::REGULAR) continue;

        auto cell = fit->first;
        int opposite = fit->second;

        std::array<int, 3> face;
        int idx = 0;
        for (int i = 0; i < 4; ++i) {
            if (i != opposite)
                face[idx++] = get_or_insert(cell->vertex(i));
        }
        faces.push_back(face);
    }

    // Volume tetrahedra: INTERIOR cells
    for (auto cit = A.finite_cells_begin(); cit != A.finite_cells_end(); ++cit) {
        if (A.classify(cit) != Fixed_alpha_shape::INTERIOR) continue;

        std::array<int, 4> tet;
        for (int i = 0; i < 4; i++)
            tet[i] = get_or_insert(cit->vertex(i));
        tets.push_back(tet);
    }

    size_t nv = vertices.size();
    size_t nf = faces.size();
    size_t nt = tets.size();

    py::array_t<float>   verts_arr({nv, size_t(3)});
    py::array_t<int32_t> faces_arr({nf, size_t(3)});
    py::array_t<int32_t> tets_arr ({nt, size_t(4)});

    auto vb = verts_arr.mutable_unchecked<2>();
    auto fb = faces_arr.mutable_unchecked<2>();
    auto tb = tets_arr.mutable_unchecked<2>();

    for (size_t i = 0; i < nv; i++) {
        auto p = vertices[i]->point().point();
        vb(i, 0) = static_cast<float>(p.x());
        vb(i, 1) = static_cast<float>(p.y());
        vb(i, 2) = static_cast<float>(p.z());
    }
    for (size_t i = 0; i < nf; i++) {
        fb(i, 0) = faces[i][0]; fb(i, 1) = faces[i][1]; fb(i, 2) = faces[i][2];
    }
    for (size_t i = 0; i < nt; i++) {
        tb(i, 0) = tets[i][0]; tb(i, 1) = tets[i][1];
        tb(i, 2) = tets[i][2]; tb(i, 3) = tets[i][3];
    }

    return std::make_tuple(verts_arr, faces_arr, tets_arr);
}

std::tuple<py::array_t<float>, py::array_t<int32_t>, py::array_t<int32_t>>
compute_alpha_complex_from_atoms(
    py::array_t<float> positions,
    py::array_t<float> radii,
    float alpha,
    float probe_radius = 1.4f
) {
    auto pos_buf = positions.unchecked<2>();
    auto rad_buf = radii.unchecked<1>();

    if (pos_buf.shape(1) != 3)
        throw std::invalid_argument("positions must have shape (N, 3)");
    if (pos_buf.shape(0) != rad_buf.shape(0))
        throw std::invalid_argument("positions and radii must have same number of atoms");

    size_t n_atoms = pos_buf.shape(0);

    std::vector<Weighted_point> weighted_points;
    weighted_points.reserve(n_atoms);
    for (size_t i = 0; i < n_atoms; i++) {
        Point_3 p(pos_buf(i, 0), pos_buf(i, 1), pos_buf(i, 2));
        double r = rad_buf(i) + probe_radius;
        weighted_points.emplace_back(p, r * r);
    }

    Triangulation T(weighted_points.begin(), weighted_points.end());
    Fixed_alpha_shape A(T, alpha);

    return extract_alpha_volume(A);
}

PYBIND11_MODULE(cgal_alpha, m) {
    m.doc() = "Alpha Complex Python bindings using CGAL Fixed_alpha_shape_3";
    m.attr("__version__") = "1.3.0";

    m.def("compute_alpha_complex_from_atoms", &compute_alpha_complex_from_atoms,
          py::arg("positions"),
          py::arg("radii"),
          py::arg("alpha"),
          py::arg("probe_radius") = 1.4f,
          R"doc(
Compute 3D alpha complex from atom positions and radii.

Args:
    positions:    (N, 3) float32 array of atom centers
    radii:        (N,)   float32 array of atom VDW radii
    alpha:        alpha value (e.g. 0.0)
    probe_radius: probe radius added to atomic radii (default 1.4 for water)

Returns:
    (vertices, faces, tetrahedra)
    vertices:     (V, 3) float32 - 3D positions
    faces:        (F, 3) int32   - boundary triangles (REGULAR facets)
    tetrahedra:   (T, 4) int32   - interior volume cells (INTERIOR cells)
)doc");
}