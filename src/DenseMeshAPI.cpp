// Copyright Gradientspace Corp. All Rights Reserved.

#include "DenseMeshAPI.h"

using namespace GS;


void DenseMeshBuilder::ResetMesh()
{
	Vertices.clear();
	Triangles.clear();
	TriangleGroups.clear();
	TriangleMaterialIDs.clear();
	Colors.clear();
	Normals.clear();
	UVs.clear();
	TriangleColorIndices.clear();
	TriangleNormalIndices.clear();
	TriangleUVIndices.clear();
	NextGroupID = 0;
}

int DenseMeshBuilder::AppendVertex(const Vector3d& Position)
{
	int Index = (int)Vertices.size();
	Vertices.push_back(Position);
	return Index;
}

int DenseMeshBuilder::GetVertexCount() const
{
	return (int)Vertices.size();
}

int DenseMeshBuilder::AllocateGroupID()
{
	return NextGroupID++;
}

int DenseMeshBuilder::AppendTriangle(const Index3i& Triangle, int GroupID)
{
	int Index = (int)Triangles.size();
	Triangles.push_back(Triangle);
	TriangleGroups.push_back(GroupID);
	TriangleMaterialIDs.push_back(0);
	TriangleColorIndices.push_back(Index3i(-1));
	TriangleNormalIndices.push_back(Index3i(-1));
	TriangleUVIndices.push_back(Index3i(-1));
	return Index;
}

int DenseMeshBuilder::GetTriangleCount() const
{
	return (int)Triangles.size();
}

void DenseMeshBuilder::SetMaterialID(int TriangleID, int MaterialID)
{
	TriangleMaterialIDs[TriangleID] = MaterialID;
}

int DenseMeshBuilder::AppendColor(const Vector4f& Color, bool bIsLinearColor)
{
	int Index = (int)Colors.size();
	Colors.push_back(Color);
	return Index;
}

void DenseMeshBuilder::SetTriangleColors(int TriangleID, const Index3i& Indices)
{
	TriangleColorIndices[TriangleID] = Indices;
}

int DenseMeshBuilder::AppendNormal(const Vector3f& Normal)
{
	int Index = (int)Normals.size();
	Normals.push_back(Normal);
	return Index;
}

void DenseMeshBuilder::SetTriangleNormals(int TriangleID, const Index3i& Indices)
{
	TriangleNormalIndices[TriangleID] = Indices;
}

int DenseMeshBuilder::AppendUV(const Vector2f& UV)
{
	int Index = (int)UVs.size();
	UVs.push_back(UV);
	return Index;
}

void DenseMeshBuilder::SetTriangleUVs(int TriangleID, const Index3i& Indices)
{
	TriangleUVIndices[TriangleID] = Indices;
}


// ----- DenseMeshCollector -----

void DenseMeshCollector::AppendMesh(const IMeshBuilder* Builder)
{
	const DenseMeshBuilder* Mesh = static_cast<const DenseMeshBuilder*>(Builder);

	int VertexOffset = (int)AccumulatedMesh.Vertices.size();
	int ColorOffset = (int)AccumulatedMesh.Colors.size();
	int NormalOffset = (int)AccumulatedMesh.Normals.size();
	int UVOffset = (int)AccumulatedMesh.UVs.size();
	int GroupOffset = AccumulatedMesh.NextGroupID;

	// append vertices
	AccumulatedMesh.Vertices.insert(AccumulatedMesh.Vertices.end(),
		Mesh->Vertices.begin(), Mesh->Vertices.end());

	// append attribute arrays
	AccumulatedMesh.Colors.insert(AccumulatedMesh.Colors.end(),
		Mesh->Colors.begin(), Mesh->Colors.end());
	AccumulatedMesh.Normals.insert(AccumulatedMesh.Normals.end(),
		Mesh->Normals.begin(), Mesh->Normals.end());
	AccumulatedMesh.UVs.insert(AccumulatedMesh.UVs.end(),
		Mesh->UVs.begin(), Mesh->UVs.end());

	// append triangles with offset vertex indices
	for (int i = 0; i < (int)Mesh->Triangles.size(); i++)
	{
		const Index3i& Tri = Mesh->Triangles[i];
		AccumulatedMesh.Triangles.push_back(
			Index3i(Tri.A + VertexOffset, Tri.B + VertexOffset, Tri.C + VertexOffset));

		AccumulatedMesh.TriangleGroups.push_back(Mesh->TriangleGroups[i] + GroupOffset);
		AccumulatedMesh.TriangleMaterialIDs.push_back(Mesh->TriangleMaterialIDs[i]);

		// offset color indices, preserving -1 sentinel
		const Index3i& CI = Mesh->TriangleColorIndices[i];
		AccumulatedMesh.TriangleColorIndices.push_back(Index3i(
			CI.A >= 0 ? CI.A + ColorOffset : -1,
			CI.B >= 0 ? CI.B + ColorOffset : -1,
			CI.C >= 0 ? CI.C + ColorOffset : -1));

		// offset normal indices, preserving -1 sentinel
		const Index3i& NI = Mesh->TriangleNormalIndices[i];
		AccumulatedMesh.TriangleNormalIndices.push_back(Index3i(
			NI.A >= 0 ? NI.A + NormalOffset : -1,
			NI.B >= 0 ? NI.B + NormalOffset : -1,
			NI.C >= 0 ? NI.C + NormalOffset : -1));

		// offset UV indices, preserving -1 sentinel
		const Index3i& UI = Mesh->TriangleUVIndices[i];
		AccumulatedMesh.TriangleUVIndices.push_back(Index3i(
			UI.A >= 0 ? UI.A + UVOffset : -1,
			UI.B >= 0 ? UI.B + UVOffset : -1,
			UI.C >= 0 ? UI.C + UVOffset : -1));
	}

	AccumulatedMesh.NextGroupID = GroupOffset + Mesh->NextGroupID;
}

AxisBox3d DenseMeshCollector::GetBounds() const
{
	if (AccumulatedMesh.Vertices.empty())
		return AxisBox3d::Empty();

	AxisBox3d Bounds(AccumulatedMesh.Vertices[0], AccumulatedMesh.Vertices[0]);
	for (size_t i = 1; i < AccumulatedMesh.Vertices.size(); i++)
	{
		Bounds.Contain(AccumulatedMesh.Vertices[i]);
	}
	return Bounds;
}


// ----- DenseMeshBuilder::ToDenseMesh -----

DenseMesh DenseMeshBuilder::ToDenseMesh() const
{
	int NumVerts = (int)Vertices.size();
	int NumTris = (int)Triangles.size();

	DenseMesh Mesh;
	Mesh.Resize(NumVerts, NumTris);

	for (int i = 0; i < NumVerts; i++)
	{
		Mesh.SetPosition(i, Vertices[i]);
	}

	for (int i = 0; i < NumTris; i++)
	{
		Mesh.SetTriangle(i, Triangles[i]);
		Mesh.SetTriGroup(i, TriangleGroups[i]);
		Mesh.SetTriMaterialIndex(i, TriangleMaterialIDs[i]);

		// normals: look up indexed attributes into per-triangle tuple
		const Index3i& NI = TriangleNormalIndices[i];
		if (NI.A >= 0)
		{
			Mesh.SetTriVtxNormals(i, TriVtxNormals(
				Normals[NI.A], Normals[NI.B], Normals[NI.C]));
		}

		// UVs
		const Index3i& UI = TriangleUVIndices[i];
		if (UI.A >= 0)
		{
			Mesh.SetTriVtxUVs(i, TriVtxUVs(
				UVs[UI.A], UVs[UI.B], UVs[UI.C]));
		}

		// colors: convert Vector4f -> Color4b
		const Index3i& CI = TriangleColorIndices[i];
		if (CI.A >= 0)
		{
			Mesh.SetTriVtxColors(i, TriVtxColors(
				Color4b(Colors[CI.A]), Color4b(Colors[CI.B]), Color4b(Colors[CI.C])));
		}
	}

	return Mesh;
}


// ----- DenseMeshBuilderFactory -----

IMeshBuilder* DenseMeshBuilderFactory::Allocate()
{
	return new DenseMeshBuilder();
}
