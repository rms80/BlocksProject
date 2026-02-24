// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include <algorithm>
#include <cmath>

#include "Math/GSIntVector3.h"
#include "Math/GSIndex3.h"
#include "src/DenseMeshAPI.h"

namespace GS
{

template<int NumVertices>
struct MeshCode
{
	// 3 ints per vertex (x,y,z)
	static constexpr int DataSize = NumVertices * 3;
	int Data[DataSize];

	bool operator==(const MeshCode& Other) const
	{
		for (int i = 0; i < DataSize; i++)
		{
			if (Data[i] != Other.Data[i]) return false;
		}
		return true;
	}

	bool operator<(const MeshCode& Other) const
	{
		for (int i = 0; i < DataSize; i++)
		{
			if (Data[i] != Other.Data[i]) return Data[i] < Other.Data[i];
		}
		return false;
	}
};

using RampMeshCode = MeshCode<6>;


template<int NumVertices>
MeshCode<NumVertices> ComputeMeshCode(const DenseMeshBuilder& Builder)
{
	MeshCode<NumVertices> Code;

	// 1. Quantize vertices: round each component to int
	std::vector<Vector3i> QuantizedVertices(Builder.Vertices.size());
	for (size_t i = 0; i < Builder.Vertices.size(); i++)
	{
		const Vector3d& V = Builder.Vertices[i];
		QuantizedVertices[i] = Vector3i(
			(int)std::round(V.X),
			(int)std::round(V.Y),
			(int)std::round(V.Z)
		);
	}

	// 2. Build sorted index order for vertices
	std::vector<int> SortedIndices(QuantizedVertices.size());
	for (size_t i = 0; i < SortedIndices.size(); i++)
	{
		SortedIndices[i] = (int)i;
	}
	std::sort(SortedIndices.begin(), SortedIndices.end(), [&](int A, int B) {
		return QuantizedVertices[A] < QuantizedVertices[B];
	});

	// 3. Pack into MeshCode::Data
	int Idx = 0;
	for (size_t i = 0; i < SortedIndices.size(); i++)
	{
		const Vector3i& V = QuantizedVertices[SortedIndices[i]];
		Code.Data[Idx++] = V.X;
		Code.Data[Idx++] = V.Y;
		Code.Data[Idx++] = V.Z;
	}

	return Code;
}


} // end namespace GS
