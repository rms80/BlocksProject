// Copyright Gradientspace Corp. All Rights Reserved.
#pragma once

#include "Mesh/GenericMeshAPI.h"

namespace GS
{

class BoxGenerator
{
public:
	Vector3d Center = Vector3d::Zero();
	Vector3d Dimensions = Vector3d(1.0, 1.0, 1.0);

	void Generate(IMeshBuilder& Builder) const;
};

} // end namespace GS
