// Copyright Gradientspace Corp. All Rights Reserved.
#include "BoxGenerator.h"

using namespace GS;

void BoxGenerator::Generate(IMeshBuilder& Builder) const
{
	Vector3d Half = Dimensions * 0.5;
	Vector3d Min = Center - Half;
	Vector3d Max = Center + Half;

	int GroupID = Builder.AllocateGroupID();

	// 8 vertices
	int V0 = Builder.AppendVertex(Vector3d(Min.X, Min.Y, Min.Z));
	int V1 = Builder.AppendVertex(Vector3d(Max.X, Min.Y, Min.Z));
	int V2 = Builder.AppendVertex(Vector3d(Max.X, Max.Y, Min.Z));
	int V3 = Builder.AppendVertex(Vector3d(Min.X, Max.Y, Min.Z));
	int V4 = Builder.AppendVertex(Vector3d(Min.X, Min.Y, Max.Z));
	int V5 = Builder.AppendVertex(Vector3d(Max.X, Min.Y, Max.Z));
	int V6 = Builder.AppendVertex(Vector3d(Max.X, Max.Y, Max.Z));
	int V7 = Builder.AppendVertex(Vector3d(Min.X, Max.Y, Max.Z));

	// 8 normals (one per vertex)
	for (int i = 0; i < 8; ++i)
		Builder.AppendNormal(Vector3f(0, 0, 0));

	// 12 triangles (2 per face, outward-facing)
	// -Z face (V0,V2,V1), (V0,V3,V2)
	Builder.AppendTriangle(Index3i(V0, V2, V1), GroupID);
	Builder.AppendTriangle(Index3i(V0, V3, V2), GroupID);
	// +Z face (V4,V5,V6), (V4,V6,V7)
	Builder.AppendTriangle(Index3i(V4, V5, V6), GroupID);
	Builder.AppendTriangle(Index3i(V4, V6, V7), GroupID);
	// -Y face (V0,V1,V5), (V0,V5,V4)
	Builder.AppendTriangle(Index3i(V0, V1, V5), GroupID);
	Builder.AppendTriangle(Index3i(V0, V5, V4), GroupID);
	// +Y face (V2,V3,V7), (V2,V7,V6)
	Builder.AppendTriangle(Index3i(V2, V3, V7), GroupID);
	Builder.AppendTriangle(Index3i(V2, V7, V6), GroupID);
	// -X face (V0,V4,V7), (V0,V7,V3)
	Builder.AppendTriangle(Index3i(V0, V4, V7), GroupID);
	Builder.AppendTriangle(Index3i(V0, V7, V3), GroupID);
	// +X face (V1,V2,V6), (V1,V6,V5)
	Builder.AppendTriangle(Index3i(V1, V2, V6), GroupID);
	Builder.AppendTriangle(Index3i(V1, V6, V5), GroupID);
}
