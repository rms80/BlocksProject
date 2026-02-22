// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include <vector>
#include "Mesh/GenericMeshAPI.h"
#include "Mesh/DenseMesh.h"

namespace GS
{

class DenseMeshBuilder : public IMeshBuilder
{
public:
	std::vector<Vector3d> Vertices;
	std::vector<Index3i> Triangles;

	std::vector<int> TriangleGroups;
	std::vector<int> TriangleMaterialIDs;

	std::vector<Vector4f> Colors;
	std::vector<Vector3f> Normals;
	std::vector<Vector2f> UVs;

	std::vector<Index3i> TriangleColorIndices;
	std::vector<Index3i> TriangleNormalIndices;
	std::vector<Index3i> TriangleUVIndices;

	int NextGroupID = 0;

	// IMeshBuilder interface
	void ResetMesh() override;

	int AppendVertex(const Vector3d& Position) override;
	int GetVertexCount() const override;

	int AllocateGroupID() override;
	int AppendTriangle(const Index3i& Triangle, int GroupID) override;
	int GetTriangleCount() const override;

	void SetMaterialID(int TriangleID, int MaterialID) override;

	int AppendColor(const Vector4f& Color, bool bIsLinearColor) override;
	void SetTriangleColors(int TriangleID, const Index3i& TriColorIndices) override;

	int AppendNormal(const Vector3f& Normal) override;
	void SetTriangleNormals(int TriangleID, const Index3i& TriNormalIndices) override;

	int AppendUV(const Vector2f& UV) override;
	void SetTriangleUVs(int TriangleID, const Index3i& TriUVIndices) override;

	// Helpers
	bool HasColors() const { return !Colors.empty(); }
	bool HasNormals() const { return !Normals.empty(); }
	bool HasUVs() const { return !UVs.empty(); }

	// Conversion
	DenseMesh ToDenseMesh() const;
};


class DenseMeshCollector : public IMeshCollector
{
public:
	DenseMeshBuilder AccumulatedMesh;

	void AppendMesh(const IMeshBuilder* Builder) override;
	AxisBox3d GetBounds() const override;
};


class DenseMeshBuilderFactory : public IMeshBuilderFactory
{
public:
	IMeshBuilder* Allocate() override;
};

}
