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

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, const Matrix4x4 &toWorld, std::vector<hitable *> &lights)
{
    std::vector<std::shared_ptr<triangle_mesh>> meshes = mesh_loader::load_obj(file);
    std::vector<std::shared_ptr<hitable>> triangles;
    for (auto &mesh : meshes)
    {
        for (int i = 0; i < mesh->nTriangles; ++i)
        {
            // Transform vertices and normals to global coordinates using toWorld matrix
            Vector4f transformed;
            for (int j = 0; j < mesh->nVertices; ++j)
            {
                transformed = toWorld * Vector4f(mesh->vertices[j][0], mesh->vertices[j][1], mesh->vertices[j][2], 1.0f);
                mesh->vertices[j] = Vector3f(transformed[0], transformed[1], transformed[2]);

                transformed = toWorld * Vector4f(mesh->normals[j][0], mesh->normals[j][1], mesh->normals[j][2], 1.0f);
                mesh->normals[j] = Vector3f(transformed[0], transformed[1], transformed[2]);
            }

            triangles.push_back(std::make_shared<triangle>(mesh, i));
            if (dynamic_cast<diffuse_light *>(mesh->mat.get()))
                lights.push_back(triangles.back().get());
        }
    }

    return triangles;
}
