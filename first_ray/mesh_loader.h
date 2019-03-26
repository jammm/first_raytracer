#ifndef MESH_LOADER_H_
#define MESH_LOADER_H_

#include "triangle.h"
#include <memory>
#include <assert.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/DefaultLogger.hpp>
#include <assimp/LogStream.hpp>

class mesh_loader
{
public:
    static std::shared_ptr<triangle_mesh> load_obj(std::string file);

    static std::unique_ptr<const aiScene> scene;
};

#endif
