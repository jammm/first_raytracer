#include "mesh_loader.h"

std::vector<std::shared_ptr<triangle_mesh>> mesh_loader::load_obj(std::string file)
{
    std::vector<std::shared_ptr<triangle_mesh>> meshes;

    Assimp::Importer importer;
    const aiScene *scene = importer.ReadFile(file, aiProcess_Triangulate | aiProcess_GenSmoothNormals | aiProcess_JoinIdenticalVertices);
    // If the import failed, report it
    assert(scene);

	for (unsigned int i = 0; i < scene->mNumMeshes; ++i)
	{
		aiMesh *mesh = scene->mMeshes[i];

		Vector3f *vertices = new Vector3f[mesh->mNumVertices];
		Vector3f *normals = new Vector3f[mesh->mNumVertices];
		Point2f *uv = new Point2f[mesh->mNumVertices];
        std::vector<int> indices;

		for (unsigned int j = 0; j < mesh->mNumVertices; ++j)
		{
			vertices[j][0] = mesh->mVertices[j].x;
			vertices[j][1] = mesh->mVertices[j].y;
			vertices[j][2] = mesh->mVertices[j].z;

			normals[j][0] = mesh->mNormals[j].x;
			normals[j][1] = mesh->mNormals[j].y;
			normals[j][2] = mesh->mNormals[j].z;

            // Some meshes don't have texture coords
            if (mesh->HasTextureCoords(j))
            {
                uv[j].x = mesh->mTextureCoords[0][j].x;
                uv[j].y = mesh->mTextureCoords[0][j].y;
            }
		}
		//Store indices
		for (unsigned int k = 0; k < mesh->mNumFaces; k++)
		{
			const aiFace& face = mesh->mFaces[k];
			assert(face.mNumIndices == 3);
			indices.push_back(face.mIndices[0]);
			indices.push_back(face.mIndices[1]);
			indices.push_back(face.mIndices[2]);
		}

        //Push triangle mesh into vector
        meshes.push_back(std::make_shared<triangle_mesh>(mesh->mNumFaces, mesh->mNumVertices, vertices, indices.data(), normals, uv, true));
	}

    return meshes;
}
