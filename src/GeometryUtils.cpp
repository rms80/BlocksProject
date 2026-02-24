// Copyright Gradientspace Corp. All Rights Reserved.
#include "GeometryUtils.h"
#include "Math/GSMath.h"
#include "MeshIO/OBJWriter.h"
#include "Core/TextIO.h"

#include <cmath>
#include <iostream>

using namespace GS;


double GS::TriDistanceSqr(const Vector3d& V0, const Vector3d& V1, const Vector3d& V2, const Vector3d& Point)
{
	Vector3d diff = V0 - Point;
	Vector3d edge0 = V1 - V0;
	Vector3d edge1 = V2 - V0;
	double a00 = edge0.SquaredLength();
	double a01 = edge0.Dot(edge1);
	double a11 = edge1.SquaredLength();
	double b0 = diff.Dot(edge0);
	double b1 = diff.Dot(edge1);
	double c = diff.SquaredLength();
	double det = std::abs(a00 * a11 - a01 * a01);
	double s = a01 * b1 - a11 * b0;
	double t = a01 * b0 - a00 * b1;
	double sqrDistance;

	if (s + t <= det) {
		if (s < 0) {
			if (t < 0) {
				// region 4
				if (b0 < 0) {
					t = 0;
					if (-b0 >= a00) {
						s = 1;
						sqrDistance = a00 + 2.0 * b0 + c;
					} else {
						s = -b0 / a00;
						sqrDistance = b0 * s + c;
					}
				} else {
					s = 0;
					if (b1 >= 0) {
						t = 0;
						sqrDistance = c;
					} else if (-b1 >= a11) {
						t = 1;
						sqrDistance = a11 + 2.0 * b1 + c;
					} else {
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			} else {
				// region 3
				s = 0;
				if (b1 >= 0) {
					t = 0;
					sqrDistance = c;
				} else if (-b1 >= a11) {
					t = 1;
					sqrDistance = a11 + 2.0 * b1 + c;
				} else {
					t = -b1 / a11;
					sqrDistance = b1 * t + c;
				}
			}
		} else if (t < 0) {
			// region 5
			t = 0;
			if (b0 >= 0) {
				s = 0;
				sqrDistance = c;
			} else if (-b0 >= a00) {
				s = 1;
				sqrDistance = a00 + 2.0 * b0 + c;
			} else {
				s = -b0 / a00;
				sqrDistance = b0 * s + c;
			}
		} else {
			// region 0 (interior)
			double invDet = 1.0 / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
				t * (a01 * s + a11 * t + 2.0 * b1) + c;
		}
	} else {
		double tmp0, tmp1, numer, denom;
		if (s < 0) {
			// region 2
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if (tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a00 - 2.0 * a01 + a11;
				if (numer >= denom) {
					s = 1;
					t = 0;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			} else {
				s = 0;
				if (tmp1 <= 0) {
					t = 1;
					sqrDistance = a11 + 2.0 * b1 + c;
				} else if (b1 >= 0) {
					t = 0;
					sqrDistance = c;
				} else {
					t = -b1 / a11;
					sqrDistance = b1 * t + c;
				}
			}
		} else if (t < 0) {
			// region 6
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if (tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a00 - 2.0 * a01 + a11;
				if (numer >= denom) {
					t = 1;
					s = 0;
					sqrDistance = a11 + 2.0 * b1 + c;
				} else {
					t = numer / denom;
					s = 1 - t;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			} else {
				t = 0;
				if (tmp1 <= 0) {
					s = 1;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else if (b0 >= 0) {
					s = 0;
					sqrDistance = c;
				} else {
					s = -b0 / a00;
					sqrDistance = b0 * s + c;
				}
			}
		} else {
			// region 1
			numer = a11 + b1 - a01 - b0;
			if (numer <= 0) {
				s = 0;
				t = 1;
				sqrDistance = a11 + 2.0 * b1 + c;
			} else {
				denom = a00 - 2.0 * a01 + a11;
				if (numer >= denom) {
					s = 1;
					t = 0;
					sqrDistance = a00 + 2.0 * b0 + c;
				} else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
						t * (a01 * s + a11 * t + 2.0 * b1) + c;
				}
			}
		}
	}

	return (sqrDistance < 0) ? 0 : sqrDistance;
}


double GS::TriSolidAngle(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& p)
{
	Vector3d va = a - p, vb = b - p, vc = c - p;
	double la = va.Length(), lb = vb.Length(), lc = vc.Length();
	double numerator = va.Dot(GS::Cross(vb, vc));
	double denominator = la * lb * lc + va.Dot(vb) * lc + vb.Dot(vc) * la + va.Dot(vc) * lb;
	return 2.0 * std::atan2(numerator, denominator);
}


AxisBox3d GS::ComputeBounds(const DenseMesh& Mesh)
{
	int NumVerts = Mesh.GetVertexCount();
	if (NumVerts == 0)
		return AxisBox3d::Empty();
	AxisBox3d Bounds(Mesh.GetPosition(0), Mesh.GetPosition(0));
	for (int i = 1; i < NumVerts; ++i)
		Bounds.Contain(Mesh.GetPosition(i));
	return Bounds;
}


bool GS::WriteMeshOBJ(const std::string& Path, const DenseMesh& Mesh, bool bReverseOrientation)
{
	FileTextWriter Writer = FileTextWriter::OpenFile(Path);
	if (!Writer.IsOpen()) {
		std::cerr << "Failed to open " << Path << " for writing" << std::endl;
		return false;
	}
	OBJWriter::WriteOptions Opts;
	Opts.bReverseTriOrientation = bReverseOrientation;
	Opts.bVertexColors = false;
	OBJWriter::WriteOBJ(Writer, Mesh, Opts);
	std::cout << "Wrote " << Path << std::endl;
	return true;
}


bool GS::TestBoxOverlap(const AxisBox3d& Box, FunctionRef<bool(const Vector3d&)> Predicate)
{
	Vector3d Mid = (Box.Min + Box.Max) * 0.5;

	// Test 8 corners
	for (int cz = 0; cz <= 1; ++cz)
		for (int cy = 0; cy <= 1; ++cy)
			for (int cx = 0; cx <= 1; ++cx) {
				Vector3d Corner(
					cx ? Box.Max.X : Box.Min.X,
					cy ? Box.Max.Y : Box.Min.Y,
					cz ? Box.Max.Z : Box.Min.Z);
				if (Predicate(Corner))
					return true;
			}

	// Test box center
	if (Predicate(Mid))
		return true;

	// Test 6 face centers
	if (Predicate(Vector3d(Box.Min.X, Mid.Y, Mid.Z))) return true;
	if (Predicate(Vector3d(Box.Max.X, Mid.Y, Mid.Z))) return true;
	if (Predicate(Vector3d(Mid.X, Box.Min.Y, Mid.Z))) return true;
	if (Predicate(Vector3d(Mid.X, Box.Max.Y, Mid.Z))) return true;
	if (Predicate(Vector3d(Mid.X, Mid.Y, Box.Min.Z))) return true;
	if (Predicate(Vector3d(Mid.X, Mid.Y, Box.Max.Z))) return true;

	// Test 12 edge midpoints (4 per axis)
	for (int a = 0; a <= 1; ++a)
		for (int b = 0; b <= 1; ++b) {
			double ey = a ? Box.Max.Y : Box.Min.Y;
			double ez = b ? Box.Max.Z : Box.Min.Z;
			if (Predicate(Vector3d(Mid.X, ey, ez))) return true;
		}
	for (int a = 0; a <= 1; ++a)
		for (int b = 0; b <= 1; ++b) {
			double ex = a ? Box.Max.X : Box.Min.X;
			double ez = b ? Box.Max.Z : Box.Min.Z;
			if (Predicate(Vector3d(ex, Mid.Y, ez))) return true;
		}
	for (int a = 0; a <= 1; ++a)
		for (int b = 0; b <= 1; ++b) {
			double ex = a ? Box.Max.X : Box.Min.X;
			double ey = b ? Box.Max.Y : Box.Min.Y;
			if (Predicate(Vector3d(ex, ey, Mid.Z))) return true;
		}

	return false;
}


double GS::WindingNumber(const DenseMesh& Mesh, const Vector3d& Point)
{
	static constexpr double FourPi = 4.0 * GS::RealConstants<double>::Pi();
	double Sum = 0;
	int NumTris = Mesh.GetTriangleCount();
	for (int i = 0; i < NumTris; ++i)
	{
		const Index3i& Tri = Mesh.GetTriangle(i);
		Sum += TriSolidAngle(
			Mesh.GetPosition(Tri.A), Mesh.GetPosition(Tri.B), Mesh.GetPosition(Tri.C), Point);
	}
	return Sum / FourPi;
}
