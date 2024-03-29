#include "preprocess.hpp"

int main()
{
    Preprocess preprocess;
    const Preprocess::Mesh mesh = preprocess.read_gmsh("../FEM-2D_plane_elasticity/Square.msh");
    const Preprocess::Info info = preprocess.read_info("../FEM-2D_plane_elasticity/input.yml");
    std::vector<std::vector<int>> ID = preprocess.get_ID(mesh, info);
    std::vector<std::vector<int>> IEN = preprocess.get_IEN(mesh, info);
    std::vector<std::vector<int>> LM = preprocess.get_LM(ID, IEN);
    const Preprocess::Tri_quad tri_quad = preprocess.TriQuad(info.quad_deg);
    const Preprocess::Line_quad line_quad = preprocess.Gauss(info.quad_deg, -1.0, 1.0);
    return 0;
}
