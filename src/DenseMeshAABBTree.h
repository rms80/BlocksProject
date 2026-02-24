// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Mesh/DenseMesh.h"
#include "Math/GSVector3.h"
#include "Math/GSAxisBox3.h"
#include "Math/GSIndex2.h"
#include "Core/unsafe_vector.h"

namespace GS
{

class DenseMeshAABBTree
{
public:
	using BoxType = AxisBox3d;

	struct ChildIndex {
		uint32_t LeafCount : 4;		// max 16 triangles in a leaf. If > 0, Index is a leaf-node-index
		uint32_t Index : 28;		// max 268M internal or leaf nodes
	};

	struct InteriorNode {
		ChildIndex LeftChild;
		ChildIndex RightChild;
		BoxType LeftBounds;
		BoxType RightBounds;
	};

	struct LeafTriangle {
		AxisBox3d Box;
		int32_t TriangleID;
	};

	ChildIndex RootIndex;
	BoxType RootBounds;
	unsafe_vector<InteriorNode> NodeTree;
	unsafe_vector<LeafTriangle> LeafTriLists;

	const DenseMesh* Mesh = nullptr;

	void Build(const DenseMesh& MeshIn);

	// returns triangle ID of nearest triangle, or -1 if none found
	// NearestDistSqr is set to the squared distance to the nearest triangle
	int FindNearestTriangle(const Vector3d& Point, double& NearestDistSqr, double MaxDist = 1e18) const;

	// returns the winding number at the given point using fast multipole approximation
	double FastWindingNumber(const Vector3d& Point) const;

	// returns true if the point is inside the mesh (winding number > 0.5)
	bool IsInside(const Vector3d& Point) const;

private:

	// Fast Winding Number cache (one entry per interior node)
	struct FWNInfo {
		Vector3d Center;
		double Radius;
		Vector3d Order1Vec;		// sum of area-weighted normals
	};
	unsafe_vector<FWNInfo> FWNCache;
	bool bFWNCacheBuilt = false;
	double FWNBeta = 2.0;

	void BuildFWNCache();
	void BuildFWNCacheRecursive(ChildIndex Index,
		Vector3d& CenterOut, double& RadiusOut, Vector3d& Order1VecOut, double& TotalAreaOut);

	double EvalFastWindingNumber(ChildIndex Index, const Vector3d& Point) const;
};

} // end namespace GS
