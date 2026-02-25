// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Math/GSVector3.h"
#include "Math/GSAxisBox3.h"
#include "Core/FunctionRef.h"
#include "Mesh/DenseMesh.h"

#include <string>
#include <array>

#include "DenseMeshAPI.h"

namespace GS
{

// Fan-triangulate a face polygon, removing consecutive duplicate vertex indices (from merged vertices).
// N is the max polygon size. CW winding is preserved.
template<size_t N>
void FanTriangulate(DenseMeshBuilder& Builder, int GroupID, const std::array<int, N>& Poly)
{
    std::array<int, N> P;
    size_t Count = 0;
    for (size_t i = 0; i < N; ++i) {
        if (Count == 0 || Poly[i] != P[Count - 1])
            P[Count++] = Poly[i];
    }
    if (Count > 1 && P[0] == P[Count - 1])
        Count--;
    for (size_t i = 1; i + 1 < Count; ++i)
        Builder.AppendTriangle(Index3i(P[0], P[i], P[i+1]), GroupID);
}

// Point-to-triangle squared distance
// Port of MeshQueries.TriDistanceSqr from geometry3Sharp
double TriDistanceSqr(const Vector3d& V0, const Vector3d& V1, const Vector3d& V2, const Vector3d& Point);

// Solid angle subtended by triangle (a,b,c) at point p
// Van Oosterom and Strackee (1983)
double TriSolidAngle(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& p);

// Analytic winding number of a closed mesh at a point (exact, brute-force over all triangles)
double WindingNumber(const DenseMesh& Mesh, const Vector3d& Point);

// Compute axis-aligned bounding box of a DenseMesh
AxisBox3d ComputeBounds(const DenseMesh& Mesh);

// Write a DenseMesh to an OBJ file. Returns true on success.
bool WriteMeshOBJ(const std::string& Path, const DenseMesh& Mesh, bool bReverseOrientation);

// Returns true if the given predicate is true at any corner of the bounding box
bool TestBoxOverlap(const AxisBox3d& Box, FunctionRef<bool(const Vector3d&)> Predicate);

// Check if a DenseMesh contains any degenerate triangles (any edge shorter than Epsilon).
// Returns the number of degenerate triangles found.
int CountDegenerateTriangles(const DenseMesh& Mesh, double Epsilon = 1e-8);

} // end namespace GS
