#include "triangle.h"
#include "mesh_loader.h"


std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, std::vector<hitable *> &lights)
{
    std::vector<std::shared_ptr<triangle_mesh>> meshes = mesh_loader::load_obj(file);
    std::vector<std::shared_ptr<hitable>> triangles;
    for (auto &mesh : meshes)
    {
        for (int i = 0; i < mesh->nTriangles; ++i)
        {
            triangles.push_back(std::make_shared<triangle>(mesh, i));
            if (dynamic_cast<diffuse_light *>(mesh->mat.get()))
                lights.push_back(triangles.back().get());
        }
    }

    return triangles;
}
