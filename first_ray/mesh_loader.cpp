#include "mesh_loader.h"
#include "material.h"
#include "util.h"

std::vector<std::shared_ptr<triangle_mesh>> mesh_loader::load_obj(std::string file, bool use_geometry_normals)
{
    std::vector<std::shared_ptr<triangle_mesh>> meshes;

    Assimp::Importer importer;
	const aiScene *scene = importer.ReadFile(file, (aiProcessPreset_TargetRealtime_Fast & ~aiProcess_GenNormals) | aiProcess_GenSmoothNormals);
    // If the import failed, report it
    if (!scene)
    {
        std::cout << "Failed to load scene! Error: " << importer.GetErrorString();
    }
    assert(scene);

	for (unsigned int i = 0; i < scene->mNumMeshes; ++i)
	{
		aiMesh *mesh = scene->mMeshes[i];

		Vector3f *vertices = new Vector3f[mesh->mNumVertices];
		Vector3f *normals = new Vector3f[mesh->mNumVertices];
		Vector2f *uv = new Vector2f[mesh->mNumVertices];
        std::vector<int> indices;

		for (unsigned int j = 0; j < mesh->mNumVertices; ++j)
		{
			vertices[j][0] = mesh->mVertices[j].x;
			vertices[j][1] = mesh->mVertices[j].y;
			vertices[j][2] = mesh->mVertices[j].z;

            // Some meshes don't have texture coords
            if (mesh->HasTextureCoords(0))
            {
                uv[j].x = mesh->mTextureCoords[0][j].x;
                uv[j].y = mesh->mTextureCoords[0][j].y;
            }
			if (mesh->HasNormals())
			{
				normals[j][0] = mesh->mNormals[j].x;
				normals[j][1] = mesh->mNormals[j].y;
				normals[j][2] = mesh->mNormals[j].z;
			}
		}
		//Store indices
		for (unsigned int k = 0; k < mesh->mNumFaces; ++k)
		{
			const aiFace& face = mesh->mFaces[k];
			assert(face.mNumIndices == 3);
			indices.push_back(face.mIndices[0]);
			indices.push_back(face.mIndices[1]);
			indices.push_back(face.mIndices[2]);
		}

        //std::unique_ptr<image> img = std::make_unique<image>("cube/default.png");
        //std::unique_ptr<material> mati = std::make_unique<lambertian>(new image_texture(img));

        std::unique_ptr<material> mat;
        aiString name;
        aiMaterial *matt = scene->mMaterials[mesh->mMaterialIndex];
        matt->Get(AI_MATKEY_NAME, name);
        if (name != aiString("DefaultMaterial"))
        {

            std::cout << "Loading " << name.C_Str() << std::endl;

            aiColor3D c(0.f, 0.f, 0.f);
            matt->Get(AI_MATKEY_COLOR_EMISSIVE, c);

            if (c != aiColor3D(0, 0, 0))
            {
                mat = std::make_unique<diffuse_light>(new constant_texture(Vector3f(c.r, c.g, c.b)));
            }
            else
            {
				matt->Get(AI_MATKEY_COLOR_SPECULAR, c);
				if (c != aiColor3D(0, 0, 0))
				{
					aiColor3D c_d(0.f, 0.f, 0.f);
					matt->Get(AI_MATKEY_COLOR_DIFFUSE, c_d);
					double opacity;
					matt->Get(AI_MATKEY_OPACITY, opacity);
					if (opacity < 1.0)
					{
						double ref_idx;
						matt->Get(AI_MATKEY_REFRACTI, ref_idx);
						mat = std::make_unique<dielectric>(ref_idx, new constant_texture(Vector3f(FromSrgb(c.r), FromSrgb(c.g), FromSrgb(c.b))),
							new constant_texture(Vector3f(1.0, 1.0, 1.0)));
					}
					else
					{
						double specular_exponent;
						matt->Get(AI_MATKEY_SHININESS, specular_exponent);
						mat = std::make_unique<modified_phong>(new constant_texture(Vector3f(FromSrgb(c_d.r), FromSrgb(c_d.g), FromSrgb(c_d.b))),
															   new constant_texture(Vector3f(FromSrgb(c.r), FromSrgb(c.g), FromSrgb(c.b))), specular_exponent);
					}
				}
				else
				{
					matt->Get(AI_MATKEY_COLOR_DIFFUSE, c);
					mat = std::make_unique<lambertian>(new constant_texture(Vector3f(FromSrgb(c.r), FromSrgb(c.g), FromSrgb(c.b))));
				}

            }
        }
        else
        {
            // Use placeholder material (lambertian with white reflectance)
            mat = std::make_unique<lambertian>(new constant_texture(Vector3f(0.5, 0.5, 0.5)));
            name = aiString("unknown");
        }

		if (!(mesh->HasNormals()))
		{
			memset(normals, 0, sizeof(Vector3f) * mesh->mNumVertices);
			// Generate well-behaved normals
			// Based on "Computing Vertex Normals from Polygonal Facets" - Grit Thuermer and Charles A. Wuethrich,
			// JGT 1998, Vol 3
			for (int j = 0; j < mesh->mNumFaces; ++j)
			{
				Vector3f n(0.0, 0.0, 0.0);
				const int* V = &indices[3 * j];
				for (int i = 0; i < 3; ++i)
				{
					const Vector3f& v0 = vertices[V[i]];
					const Vector3f& v1 = vertices[V[(i + 1) % 3]];
					const Vector3f& v2 = vertices[V[(i + 2) % 3]];
					const Vector3f edge1 = v1 - v0;
					const Vector3f edge2 = v2 - v0;

					if (i == 0)
					{
						n = cross(edge1, edge2);
						double length = n.length();
						if (length == 0) break;
						n /= length;
					}
					double angle = unit_angle(unit_vector(edge1), unit_vector(edge2));
					normals[V[i]] += n * angle;
				}

			}
			for (unsigned int j = 0; j < mesh->mNumVertices; ++j)
			{
				Vector3f &n = normals[j];
				double length = n.length();
				if (length != 0)
					n /= length;
				else
				{
					n = Vector3f(1, 0, 0);
				}
			}
		}

        meshes.push_back(std::make_shared<triangle_mesh>(mesh->mNumFaces, mesh->mNumVertices, vertices, indices.data(), normals, uv,
							std::move(mat), use_geometry_normals, name.C_Str(), true));
	}

    return meshes;
}
