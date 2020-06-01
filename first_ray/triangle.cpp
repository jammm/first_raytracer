#ifdef USE_SSE
#include "triangle_sse.hpp"
#else
#include "triangle.h"
#endif
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

std::vector<std::shared_ptr<hitable>> create_triangle_mesh(const std::string &file, const Matrix4x4 &toWorld, material *bsdf, std::vector<hitable *> &lights)
{
    std::vector<std::shared_ptr<triangle_mesh>> meshes = mesh_loader::load_obj(file);
    std::vector<std::shared_ptr<hitable>> triangles;

    Matrix4x4 toWorld_inv_transform;
    if (!toWorld.invert(toWorld_inv_transform))
        throw std::runtime_error("Cannot compute inverse. Matrix isn't singular.");
    for (auto &mesh : meshes)
    {
        // Replace default material
        if (bsdf)
            mesh->mat.reset(bsdf);
        // Transform vertices and normals to global coordinates using toWorld matrix
        Vector3f transformed;
        for (int j = 0; j < mesh->nVertices; ++j)
        {
            transformed = toWorld * mesh->vertices[j];
            mesh->vertices[j] = transformed;
            transformed = toWorld_inv_transform.normal_transform(mesh->normals[j]);
            mesh->normals[j] = transformed;
        }

        // Then push the triangles
        for (int i = 0; i < mesh->nTriangles; ++i)
        {
            triangles.push_back(std::make_shared<triangle>(mesh, i));
            if (dynamic_cast<diffuse_light *>(mesh->mat.get()))
                lights.push_back(triangles.back().get());
        }
    }

    return triangles;
}
