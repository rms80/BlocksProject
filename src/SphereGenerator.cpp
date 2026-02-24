// Copyright Gradientspace Corp. All Rights Reserved.
#include "SphereGenerator.h"
#include "Math/GSMath.h"

#include <cmath>

using namespace GS;

void SphereGenerator::Generate(IMeshBuilder& Builder) const
{
	constexpr double Pi = RealConstants<double>::Pi();
	constexpr double TwoPi = 2.0 * RealConstants<double>::Pi();

	int GroupID = Builder.AllocateGroupID();

	// vertex 0: north pole
	int NorthPole = Builder.AppendVertex(Center + Vector3d(0, Radius, 0));
	Builder.AppendNormal(Vector3f(0, 1, 0));

	// interior ring vertices (Stacks-1 rings of Slices vertices each)
	int FirstRingVertex = NorthPole + 1;
	for (int stack = 1; stack < Stacks; ++stack)
	{
		double phi = Pi * (double)stack / (double)Stacks;
		double sinPhi = std::sin(phi);
		double cosPhi = std::cos(phi);
		for (int slice = 0; slice < Slices; ++slice)
		{
			double theta = TwoPi * (double)slice / (double)Slices;
			double x = sinPhi * std::cos(theta);
			double y = cosPhi;
			double z = sinPhi * std::sin(theta);
			Builder.AppendVertex(Center + Vector3d(x * Radius, y * Radius, z * Radius));
			Builder.AppendNormal(Vector3f((float)x, (float)y, (float)z));
		}
	}

	// vertex last: south pole
	int SouthPole = Builder.AppendVertex(Center + Vector3d(0, -Radius, 0));
	Builder.AppendNormal(Vector3f(0, -1, 0));

	// north pole cap triangles
	for (int slice = 0; slice < Slices; ++slice)
	{
		int next = (slice + 1) % Slices;
		int TriID = Builder.AppendTriangle(
			Index3i(NorthPole, FirstRingVertex + next, FirstRingVertex + slice), GroupID);
		Builder.SetTriangleNormals(TriID,
			Index3i(NorthPole, FirstRingVertex + next, FirstRingVertex + slice));
	}

	// body quads (split into two triangles each)
	for (int stack = 0; stack < Stacks - 2; ++stack)
	{
		int ringStart = FirstRingVertex + stack * Slices;
		int nextRingStart = ringStart + Slices;
		for (int slice = 0; slice < Slices; ++slice)
		{
			int next = (slice + 1) % Slices;
			int a = ringStart + slice;
			int b = ringStart + next;
			int c = nextRingStart + next;
			int d = nextRingStart + slice;

			int T0 = Builder.AppendTriangle(Index3i(a, b, c), GroupID);
			Builder.SetTriangleNormals(T0, Index3i(a, b, c));

			int T1 = Builder.AppendTriangle(Index3i(a, c, d), GroupID);
			Builder.SetTriangleNormals(T1, Index3i(a, c, d));
		}
	}

	// south pole cap triangles
	int lastRingStart = FirstRingVertex + (Stacks - 2) * Slices;
	for (int slice = 0; slice < Slices; ++slice)
	{
		int next = (slice + 1) % Slices;
		int TriID = Builder.AppendTriangle(
			Index3i(SouthPole, lastRingStart + slice, lastRingStart + next), GroupID);
		Builder.SetTriangleNormals(TriID,
			Index3i(SouthPole, lastRingStart + slice, lastRingStart + next));
	}
}
