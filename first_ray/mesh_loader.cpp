#include "mesh_loader.h"

std::shared_ptr<triangle_mesh> mesh_loader::load_obj(std::string file)
{
    std::vector<int> indices;
    Vector3f *vertices;
    Vector3f *normals;
    Point2f *uv;

    Assimp::Importer importer;
    const aiScene *scene = importer.ReadFile("cube/cube.obj", aiProcessPreset_TargetRealtime_Quality);
    // If the import failed, report it
    if (!scene)
    {
        fprintf(stderr, "ASSIMP ERROR: cannot load .obj file\n");
        return NULL;
    }

    // TODO: iterate through meshes
    aiMesh *mesh = scene->mMeshes[0];

    vertices = new Vector3f[mesh->mNumVertices];
    normals = new Vector3f[mesh->mNumVertices];
    uv = new Point2f[mesh->mNumVertices];
    for (unsigned int i = 0; i < mesh->mNumVertices; i++)
    {
        vertices[i][0] = mesh->mVertices[i].x;
        vertices[i][1] = mesh->mVertices[i].y;
        vertices[i][2] = mesh->mVertices[i].z;

        normals[i][0] = mesh->mNormals[i].x;
        normals[i][1] = mesh->mNormals[i].y;
        normals[i][2] = mesh->mNormals[i].z;

        uv[i].x = mesh->mTextureCoords[0][i].x;
        uv[i].y = mesh->mTextureCoords[0][i].y;
    }
    //Store indices
    for (unsigned int i = 0; i < mesh->mNumFaces; i++)
    {
        const aiFace& face = mesh->mFaces[i];
        assert(face.mNumIndices == 3);
        indices.push_back(face.mIndices[0]);
        indices.push_back(face.mIndices[1]);
        indices.push_back(face.mIndices[2]);
    }

    return std::make_shared<triangle_mesh>(mesh->mNumFaces, mesh->mNumVertices, vertices, indices.data(), normals, uv, true);
}