// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Mesh/GenericMeshAPI.h"

namespace GS
{

class SphereGenerator
{
public:
	Vector3d Center = Vector3d::Zero();
	double Radius = 1.0;
	int Slices = 16;
	int Stacks = 12;

	void Generate(IMeshBuilder& Builder) const;
};

} // end namespace GS
