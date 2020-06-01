#ifndef MESH_LOADER_H_
#define MESH_LOADER_H_

#ifdef USE_SSE
#include "triangle_sse.hpp"
#else
#include "triangle.h"
#endif
#include <memory>
#include <vector>
#include <assert.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/DefaultLogger.hpp>
#include <assimp/LogStream.hpp>

class mesh_loader
{
public:
    static std::vector<std::shared_ptr<triangle_mesh>> load_obj(std::string file);
};

#endif
