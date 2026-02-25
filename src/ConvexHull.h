// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Math/GSVector3.h"
#include "Mesh/DenseMesh.h"

#include <array>
#include <cmath>

namespace GS
{

// Convex polyhedron defined by up to MaxPlanes half-space planes.
// A point is inside if it is on the negative side of all planes.
// Template parameter MaxPlanes should be set to the known face count
// of the cell type (e.g. 6 for box, 7 for cut corner, 10 for cylinder).
template<int MaxPlanes>
class ConvexHull
{
public:
	struct Plane
	{
		Vector3d Normal;
		double Offset;		// dot(Normal, PointOnPlane)
	};

	int NumPlanes = 0;
	std::array<Plane, MaxPlanes> Planes;

	// Build planes from a convex DenseMesh. Only adds unique planes
	// (skips triangles whose normal matches an already-seen plane).
	// Auto-detects normal orientation by testing against mesh centroid.
	void BuildFromMesh(const DenseMesh& Mesh, double NormalTolerance = 0.999)
	{
		NumPlanes = 0;
		int NumTris = Mesh.GetTriangleCount();
		for (int i = 0; i < NumTris && NumPlanes < MaxPlanes; ++i)
		{
			const Index3i& Tri = Mesh.GetTriangle(i);
			Vector3d A = Mesh.GetPosition(Tri.A);
			Vector3d B = Mesh.GetPosition(Tri.B);
			Vector3d C = Mesh.GetPosition(Tri.C);
			Vector3d CrossN = GS::Cross(B - A, C - A);
			double Len = CrossN.Length();
			if (Len < 1e-12)
				continue;
			Vector3d Normal = CrossN * (1.0 / Len);

			// Check if we already have a plane with this normal
			bool bDuplicate = false;
			for (int j = 0; j < NumPlanes; ++j)
			{
				if (Normal.Dot(Planes[j].Normal) > NormalTolerance)
				{
					bDuplicate = true;
					break;
				}
			}
			if (bDuplicate)
				continue;

			Planes[NumPlanes].Normal = Normal;
			Planes[NumPlanes].Offset = Normal.Dot(A);
			NumPlanes++;
		}

		// Auto-detect orientation: compute mesh centroid and check if it's
		// on the positive side of the first plane. If so, normals are inward-facing
		// and we need to flip them so IsInside works correctly.
		if (NumPlanes > 0)
		{
			Vector3d Centroid = Vector3d::Zero();
			int NumVerts = Mesh.GetVertexCount();
			for (int i = 0; i < NumVerts; ++i)
				Centroid = Centroid + Mesh.GetPosition(i);
			Centroid = Centroid * (1.0 / NumVerts);

			double CentroidSide = Planes[0].Normal.Dot(Centroid) - Planes[0].Offset;
			if (CentroidSide > 0)
			{
				for (int i = 0; i < NumPlanes; ++i)
				{
					Planes[i].Normal = Planes[i].Normal * -1.0;
					Planes[i].Offset = -Planes[i].Offset;
				}
			}
		}
	}

	bool IsInside(const Vector3d& P, double Epsilon = 1e-10) const
	{
		for (int i = 0; i < NumPlanes; ++i)
		{
			if (Planes[i].Normal.Dot(P) > Planes[i].Offset + Epsilon)
				return false;
		}
		return true;
	}
};

} // end namespace GS
