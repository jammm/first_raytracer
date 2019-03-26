#include "triangle.h"
#include "mesh_loader.h"

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, std::shared_ptr<material> mat)
{
    std::shared_ptr<triangle_mesh> mesh = mesh_loader::load_obj(file);
    std::vector<std::shared_ptr<hitable>> triangles;
    triangles.reserve(mesh->nTriangles);
    for (int i = 0; i < mesh->nTriangles; ++i)
    {
        triangles.push_back(std::make_shared<triangle>(mesh, i, mat));
    }

    return triangles;
}