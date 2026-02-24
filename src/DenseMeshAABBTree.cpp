// Copyright Gradientspace Corp. All Rights Reserved.
#include "DenseMeshAABBTree.h"
#include "GeometryUtils.h"
#include "Math/GSMath.h"
#include "Math/GSIndex2.h"
#include "Core/gs_debug.h"

#include <cmath>

using namespace GS;


// ---- Build helpers (mirrors AxisBoxTree2 build strategy) ----

struct SourceTriBox
{
	AxisBox3d Box;
	int32_t TriangleID;
};

struct BoxGroup3
{
	AxisBox3d Bounds;
	int32_t StartIndex = -1;
	int32_t EndIndex = -1;		// inclusive
	Index2i Children = Index2i(-1, -1);
	int32_t Count() const { return EndIndex - StartIndex + 1; }
	bool IsLeaf() const { return Children.A == -1; }
};


static int select_split_axis_longest(const AxisBox3d& Box)
{
	double X = Box.Dimension(0);
	double Y = Box.Dimension(1);
	double Z = Box.Dimension(2);
	if (X >= Y && X >= Z) return 0;
	if (Y >= Z) return 1;
	return 2;
}


static void split_group_trivial(
	unsafe_vector<SourceTriBox>& AllBoxes,
	const BoxGroup3& Group,
	BoxGroup3& LeftChild,
	BoxGroup3& RightChild)
{
	gs_debug_assert(Group.Count() > 1);

	LeftChild.StartIndex = Group.StartIndex;
	LeftChild.EndIndex = Group.StartIndex + (Group.Count() / 2);
	LeftChild.Bounds = AxisBox3d::Empty();
	for (int32_t k = LeftChild.StartIndex; k <= LeftChild.EndIndex; ++k)
		LeftChild.Bounds.Contain(AllBoxes[k].Box);
	LeftChild.Children = Index2i(-1, -1);

	RightChild.StartIndex = LeftChild.EndIndex + 1;
	RightChild.EndIndex = Group.EndIndex;
	RightChild.Bounds = AxisBox3d::Empty();
	for (int32_t k = RightChild.StartIndex; k <= RightChild.EndIndex; ++k)
		RightChild.Bounds.Contain(AllBoxes[k].Box);
	RightChild.Children = Index2i(-1, -1);
}


static bool split_group_by_midpoint(
	unsafe_vector<SourceTriBox>& AllBoxes,
	const BoxGroup3& Group,
	int32_t split_axis,
	double split_axis_origin,
	BoxGroup3& LeftChild,
	BoxGroup3& RightChild)
{
	gs_debug_assert(Group.Count() > 1);

	int32_t l = Group.StartIndex;
	int32_t r = Group.EndIndex;

	while (l < r)
	{
		double lc = AllBoxes[l].Box.Center()[split_axis];
		if (lc <= split_axis_origin) {
			l++;
			continue;
		}

		while (r > l)
		{
			double rc = AllBoxes[r].Box.Center()[split_axis];
			if (rc > split_axis_origin) {
				r--;
				continue;
			} else
				break;
		}

		if (r == l)
			break;

		GS::SwapTemp(AllBoxes[l], AllBoxes[r]);
		l++;
		r--;
	}

	if (l != r) {
		GS::SwapTemp(l, r);
	}
	else {
		Vector3d C = AllBoxes[l].Box.Center();
		if (C[split_axis] > split_axis_origin)
			l--;
	}

	if (l == Group.StartIndex || r == Group.EndIndex)
	{
		split_group_trivial(AllBoxes, Group, LeftChild, RightChild);
		return true;
	}

	LeftChild.Bounds = AllBoxes[l].Box;
	for (int k = Group.StartIndex; k < l; ++k)
		LeftChild.Bounds.Contain(AllBoxes[k].Box);
	LeftChild.StartIndex = Group.StartIndex;
	LeftChild.EndIndex = l;
	LeftChild.Children = Index2i(-1, -1);

	RightChild.Bounds = AllBoxes[Group.EndIndex].Box;
	for (int k = l + 1; k < Group.EndIndex; ++k)
		RightChild.Bounds.Contain(AllBoxes[k].Box);
	RightChild.StartIndex = l + 1;
	RightChild.EndIndex = Group.EndIndex;
	RightChild.Children = Index2i(-1, -1);

	return true;
}


// ---- Build ----

void DenseMeshAABBTree::Build(const DenseMesh& MeshIn)
{
	static_assert(sizeof(ChildIndex) == sizeof(int32_t));

	Mesh = &MeshIn;
	bFWNCacheBuilt = false;

	int NumTris = MeshIn.GetTriangleCount();
	AxisBox3d CombinedBox = AxisBox3d::Empty();
	unsafe_vector<SourceTriBox> BoxesList;
	BoxesList.reserve(NumTris);
	for (int k = 0; k < NumTris; ++k) {
		const Index3i& Tri = MeshIn.GetTriangle(k);
		const Vector3d& A = MeshIn.GetPosition(Tri.A);
		const Vector3d& B = MeshIn.GetPosition(Tri.B);
		const Vector3d& C = MeshIn.GetPosition(Tri.C);
		AxisBox3d TriBox(A, A);
		TriBox.Contain(B);
		TriBox.Contain(C);

		SourceTriBox NewBox;
		NewBox.Box = TriBox;
		NewBox.TriangleID = k;
		BoxesList.add_ref(NewBox);
		CombinedBox.Contain(TriBox);
	}
	int32_t NumBoxes = (int32_t)BoxesList.size();

	int NumChildrenInLeaf = 4;

	// trivial case
	if (NumBoxes <= NumChildrenInLeaf) {
		LeafTriLists.reserve(NumBoxes);
		for (int j = 0; j < NumBoxes; j++) {
			SourceTriBox& SrcBox = BoxesList[j];
			LeafTriLists.add_ref(LeafTriangle{ SrcBox.Box, SrcBox.TriangleID });
		}
		this->RootBounds = CombinedBox;
		this->RootIndex.Index = 0;
		this->RootIndex.LeafCount = NumBoxes;
		return;
	}

	unsafe_vector<BoxGroup3> GroupTree;
	GroupTree.reserve(NumBoxes / 2);
	BoxGroup3 RootGroup{ CombinedBox, 0, NumBoxes - 1, Index2i(-1, -1) };
	GroupTree.add_ref(RootGroup);

	unsafe_vector<int32_t> SplitJobs;
	SplitJobs.reserve(NumBoxes / 8);
	SplitJobs.add(0);

	int NextIndex = -1;
	int NumInternalNodes = 0, NumLeafNodes = 0, TotalLeafItemCount = 0;
	while (SplitJobs.pop_back(NextIndex))
	{
		BoxGroup3 CurGroup = GroupTree[NextIndex];
		if (CurGroup.Count() > NumChildrenInLeaf)
		{
			NumInternalNodes++;

			int split_axis = select_split_axis_longest(CurGroup.Bounds);
			double split_axis_origin = CurGroup.Bounds.Center()[split_axis];
			BoxGroup3 LeftChild, RightChild;
			split_group_by_midpoint(BoxesList, CurGroup, split_axis, split_axis_origin, LeftChild, RightChild);

			CurGroup.Children.A = (int)GroupTree.add_ref(LeftChild);
			CurGroup.Children.B = (int)GroupTree.add_ref(RightChild);
			GroupTree[NextIndex] = CurGroup;

			SplitJobs.add(CurGroup.Children.A);
			SplitJobs.add(CurGroup.Children.B);
		}
		else
		{
			NumLeafNodes++;
			TotalLeafItemCount += CurGroup.Count();
		}
	}
	int NumTreeNodes = (int)GroupTree.size();
	gs_debug_assert(NumLeafNodes > 0 && TotalLeafItemCount == NumBoxes);

	// reassign linear indices to internal nodes
	int LinearIndex = 0;
	for (int i = 0; i < NumTreeNodes; ++i) {
		if (!GroupTree[i].IsLeaf())
			GroupTree[i].StartIndex = LinearIndex++;
	}
	gs_debug_assert(LinearIndex == NumInternalNodes);

	// build leaf list
	LeafTriLists.reserve(TotalLeafItemCount);
	for (int i = 0; i < NumTreeNodes; ++i) {
		BoxGroup3& Group = GroupTree[i];
		if (!Group.IsLeaf())
			continue;

		int start_index = (int)LeafTriLists.size(), count = Group.Count();
		for (int j = Group.StartIndex; j <= Group.EndIndex; j++)
		{
			SourceTriBox& SrcBox = BoxesList[j];
			LeafTriLists.add_ref(LeafTriangle{ SrcBox.Box, SrcBox.TriangleID });
		}
		Group.StartIndex = start_index; Group.EndIndex = count;
	}
	gs_debug_assert(LeafTriLists.size() == (size_t)NumBoxes);

	// compact into InteriorNode tree
	this->RootIndex = ChildIndex{ 0, 0 };
	this->RootBounds = CombinedBox;
	NodeTree.reserve(NumInternalNodes);

	for (int i = 0; i < NumTreeNodes; ++i)
	{
		const BoxGroup3& Group = GroupTree[i];
		if (Group.IsLeaf())
			continue;

		InteriorNode NewInteriorNode;

		const BoxGroup3& LeftChildGroup = GroupTree[Group.Children.A];
		const BoxGroup3& RightChildGroup = GroupTree[Group.Children.B];

		NewInteriorNode.LeftBounds = LeftChildGroup.Bounds;
		NewInteriorNode.LeftChild.Index = LeftChildGroup.StartIndex;
		NewInteriorNode.LeftChild.LeafCount = LeftChildGroup.IsLeaf() ? LeftChildGroup.EndIndex : 0;

		NewInteriorNode.RightBounds = RightChildGroup.Bounds;
		NewInteriorNode.RightChild.Index = RightChildGroup.StartIndex;
		NewInteriorNode.RightChild.LeafCount = RightChildGroup.IsLeaf() ? RightChildGroup.EndIndex : 0;

		[[maybe_unused]] int NewIndex = (int)NodeTree.add_ref(NewInteriorNode);
		gs_debug_assert(NewIndex == Group.StartIndex);
	}
}


// ---- FindNearestTriangle ----

int DenseMeshAABBTree::FindNearestTriangle(const Vector3d& Point, double& NearestDistSqr, double MaxDist) const
{
	gs_debug_assert(Mesh != nullptr);

	struct StackEntry
	{
		ChildIndex Index;
		AxisBox3d Box;
	};

	unsafe_vector<StackEntry> Stack;
	Stack.reserve(32);

	NearestDistSqr = MaxDist * MaxDist;
	int NearestTriID = -1;

	double RootDistSqr = RootBounds.DistanceSquared(Point);
	if (RootDistSqr > NearestDistSqr)
		return -1;

	Stack.push_back({ RootIndex, RootBounds });
	StackEntry Next;
	while (Stack.pop_back(Next))
	{
		if (Next.Box.DistanceSquared(Point) > NearestDistSqr)
			continue;

		if (Next.Index.LeafCount > 0) {
			for (uint32_t j = 0; j < Next.Index.LeafCount; ++j) {
				const LeafTriangle& Leaf = LeafTriLists[Next.Index.Index + j];
				if (Leaf.Box.DistanceSquared(Point) < NearestDistSqr) {
					const Index3i& Tri = Mesh->GetTriangle(Leaf.TriangleID);
					double DistSqr = GS::TriDistanceSqr(
						Mesh->GetPosition(Tri.A), Mesh->GetPosition(Tri.B), Mesh->GetPosition(Tri.C), Point);
					if (DistSqr < NearestDistSqr) {
						NearestDistSqr = DistSqr;
						NearestTriID = Leaf.TriangleID;
					}
				}
			}
		}
		else
		{
			const InteriorNode& Node = NodeTree[Next.Index.Index];
			StackEntry Children[2] = {
				{ Node.LeftChild, Node.LeftBounds },
				{ Node.RightChild, Node.RightBounds }
			};
			double ChildDists[2] = {
				Node.LeftBounds.DistanceSquared(Point),
				Node.RightBounds.DistanceSquared(Point)
			};
			// push farther child first so nearer child is processed first
			if (ChildDists[0] > ChildDists[1]) {
				GS::SwapTemp(ChildDists[0], ChildDists[1]);
				GS::SwapTemp(Children[0], Children[1]);
			}
			if (ChildDists[0] < NearestDistSqr) {
				if (ChildDists[1] < NearestDistSqr)
					Stack.push_back(Children[1]);
				Stack.push_back(Children[0]);
			}
		}
	}

	return NearestTriID;
}


// ---- Fast Winding Number ----

double DenseMeshAABBTree::FastWindingNumber(const Vector3d& Point) const
{
	gs_debug_assert(Mesh != nullptr);

	if (!bFWNCacheBuilt)
		const_cast<DenseMeshAABBTree*>(this)->BuildFWNCache();

	return EvalFastWindingNumber(RootIndex, Point);
}

bool DenseMeshAABBTree::IsInside(const Vector3d& Point) const
{
	return FastWindingNumber(Point) > 0.5;
}


void DenseMeshAABBTree::BuildFWNCache()
{
	int NumNodes = (int)NodeTree.size();
	FWNCache.resize(NumNodes);

	Vector3d Center;
	double Radius, TotalArea;
	Vector3d Order1Vec;
	BuildFWNCacheRecursive(RootIndex, Center, Radius, Order1Vec, TotalArea);
	bFWNCacheBuilt = true;
}


void DenseMeshAABBTree::BuildFWNCacheRecursive(ChildIndex Index,
	Vector3d& CenterOut, double& RadiusOut, Vector3d& Order1VecOut, double& TotalAreaOut)
{
	CenterOut = Vector3d::Zero();
	Order1VecOut = Vector3d::Zero();
	TotalAreaOut = 0;
	RadiusOut = 0;

	if (Index.LeafCount > 0) {
		// leaf node: compute area-weighted centroid and order-1 vector
		for (uint32_t j = 0; j < Index.LeafCount; ++j) {
			const LeafTriangle& Leaf = LeafTriLists[Index.Index + j];
			const Index3i& Tri = Mesh->GetTriangle(Leaf.TriangleID);
			const Vector3d& A = Mesh->GetPosition(Tri.A);
			const Vector3d& B = Mesh->GetPosition(Tri.B);
			const Vector3d& C = Mesh->GetPosition(Tri.C);
			Vector3d AreaNormal = GS::Cross(B - A, C - A) * 0.5;
			double Area = AreaNormal.Length();
			Vector3d Centroid = (A + B + C) * (1.0 / 3.0);
			CenterOut = CenterOut + Centroid * Area;
			Order1VecOut = Order1VecOut + AreaNormal;
			TotalAreaOut += Area;
		}
		if (TotalAreaOut > 0)
			CenterOut = CenterOut * (1.0 / TotalAreaOut);

		// compute radius
		for (uint32_t j = 0; j < Index.LeafCount; ++j) {
			const LeafTriangle& Leaf = LeafTriLists[Index.Index + j];
			const Index3i& Tri = Mesh->GetTriangle(Leaf.TriangleID);
			RadiusOut = GS::Max(RadiusOut, CenterOut.Distance(Mesh->GetPosition(Tri.A)));
			RadiusOut = GS::Max(RadiusOut, CenterOut.Distance(Mesh->GetPosition(Tri.B)));
			RadiusOut = GS::Max(RadiusOut, CenterOut.Distance(Mesh->GetPosition(Tri.C)));
		}
		return;
	}

	// interior node: recurse into children and merge
	const InteriorNode& Node = NodeTree[Index.Index];

	Vector3d LeftCenter, RightCenter, LeftOrder1, RightOrder1;
	double LeftRadius, RightRadius, LeftArea, RightArea;

	BuildFWNCacheRecursive(Node.LeftChild, LeftCenter, LeftRadius, LeftOrder1, LeftArea);
	BuildFWNCacheRecursive(Node.RightChild, RightCenter, RightRadius, RightOrder1, RightArea);

	TotalAreaOut = LeftArea + RightArea;
	Order1VecOut = LeftOrder1 + RightOrder1;

	if (TotalAreaOut > 0)
		CenterOut = (LeftCenter * LeftArea + RightCenter * RightArea) * (1.0 / TotalAreaOut);
	else
		CenterOut = (LeftCenter + RightCenter) * 0.5;

	// radius: max distance from center to any vertex in the subtree
	// approximate using child center + child radius
	RadiusOut = GS::Max(
		CenterOut.Distance(LeftCenter) + LeftRadius,
		CenterOut.Distance(RightCenter) + RightRadius);

	// store cache for this interior node
	FWNInfo& Info = FWNCache[Index.Index];
	Info.Center = CenterOut;
	Info.Radius = RadiusOut;
	Info.Order1Vec = Order1VecOut;
}


double DenseMeshAABBTree::EvalFastWindingNumber(ChildIndex Index, const Vector3d& Point) const
{
	static constexpr double FourPi = 4.0 * GS::RealConstants<double>::Pi();

	if (Index.LeafCount > 0) {
		// leaf: sum exact solid angles
		double Sum = 0;
		for (uint32_t j = 0; j < Index.LeafCount; ++j) {
			const LeafTriangle& Leaf = LeafTriLists[Index.Index + j];
			const Index3i& Tri = Mesh->GetTriangle(Leaf.TriangleID);
			Sum += GS::TriSolidAngle(
				Mesh->GetPosition(Tri.A), Mesh->GetPosition(Tri.B), Mesh->GetPosition(Tri.C), Point);
		}
		return Sum / FourPi;
	}

	const InteriorNode& Node = NodeTree[Index.Index];
	const FWNInfo& Info = FWNCache[Index.Index];

	// check if we can use the order-1 approximation for the entire subtree
	double DistToCenter = Point.Distance(Info.Center);
	bool bContained = Node.LeftBounds.Contains(Point) || Node.RightBounds.Contains(Point);
	if (!bContained && DistToCenter > FWNBeta * Info.Radius)
	{
		// order-1 approximation
		Vector3d D = Info.Center - Point;
		double Len = D.Length();
		if (Len > 0) {
			return D.Dot(Info.Order1Vec) / (FourPi * Len * Len * Len);
		}
	}

	// recurse into children
	return EvalFastWindingNumber(Node.LeftChild, Point) + EvalFastWindingNumber(Node.RightChild, Point);
}
