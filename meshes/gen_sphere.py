#!/usr/bin/env python3
"""
Generate a UV sphere OBJ file and print it to stdout.

Sphere parameters:
  - Center: origin (0, 0, 0)
  - Radius: 1.0
  - Longitude segments (slices): 16
  - Latitude segments (stacks): 12

The mesh is a closed manifold (watertight). Pole caps use triangles,
and the body bands use quads split into two triangles each.

Vertices and normals are generated for each unique position on the sphere.
For a unit sphere centered at origin, each vertex normal equals the vertex
position (already unit length).

Vertex layout:
  - Index 1: north pole (0, 1, 0)
  - Indices 2 .. (num_lat-1)*num_lon+1: body rings from top to bottom
    Ring i (1-based, i=1..num_lat-1) has num_lon vertices at latitude angle
    pi*i/num_lat, with longitude angles 2*pi*j/num_lon for j=0..num_lon-1.
  - Last index: south pole (0, -1, 0)

Face winding is counter-clockwise when viewed from outside (standard OBJ convention).
"""

import math
import sys


def generate_sphere_obj(num_lon=16, num_lat=12, radius=1.0):
    lines = []
    lines.append("# UV Sphere: radius=1.0, lon_segments=16, lat_segments=12")
    lines.append("# Vertices: %d, Normals: %d" % (
        2 + (num_lat - 1) * num_lon,
        2 + (num_lat - 1) * num_lon,
    ))

    vertices = []
    normals = []

    # --- North pole ---
    vertices.append((0.0, radius, 0.0))
    normals.append((0.0, 1.0, 0.0))

    # --- Body rings (from top to bottom, ring index i = 1 .. num_lat-1) ---
    for i in range(1, num_lat):
        lat_angle = math.pi * i / num_lat  # 0 at north pole, pi at south pole
        sin_lat = math.sin(lat_angle)
        cos_lat = math.cos(lat_angle)
        for j in range(num_lon):
            lon_angle = 2.0 * math.pi * j / num_lon
            x = sin_lat * math.cos(lon_angle)
            y = cos_lat
            z = sin_lat * math.sin(lon_angle)
            vertices.append((x * radius, y * radius, z * radius))
            normals.append((x, y, z))

    # --- South pole ---
    vertices.append((0.0, -radius, 0.0))
    normals.append((0.0, -1.0, 0.0))

    # Emit vertices
    for v in vertices:
        lines.append("v %.6f %.6f %.6f" % v)

    # Emit normals
    for n in normals:
        lines.append("vn %.6f %.6f %.6f" % n)

    # --- Faces ---
    # OBJ indices are 1-based.
    # North pole vertex: index 1
    # Ring i (1-based) vertex j (0-based): (i-1)*num_lon + j + 2
    # South pole vertex: len(vertices)

    def ring_idx(ring, seg):
        """Return 1-based vertex index for ring `ring` (1..num_lat-1),
        segment `seg` (0..num_lon-1, wraps around)."""
        return (ring - 1) * num_lon + (seg % num_lon) + 2

    north_pole = 1
    south_pole = len(vertices)

    # North pole cap: triangles connecting north pole to the first ring.
    # Winding: north_pole, ring1[j], ring1[j+1]  (CCW from outside)
    for j in range(num_lon):
        v0 = north_pole
        v1 = ring_idx(1, j)
        v2 = ring_idx(1, j + 1)
        lines.append("f %d//%d %d//%d %d//%d" % (v0, v0, v1, v1, v2, v2))

    # Body bands: quads split into two triangles each.
    # For ring i to ring i+1 (i = 1 .. num_lat-2):
    #   Quad corners (CCW from outside):
    #     A = ring_idx(i, j)
    #     B = ring_idx(i+1, j)
    #     C = ring_idx(i+1, j+1)
    #     D = ring_idx(i, j+1)
    #   Triangle 1: A, B, D
    #   Triangle 2: B, C, D
    for i in range(1, num_lat - 1):
        for j in range(num_lon):
            a = ring_idx(i, j)
            b = ring_idx(i + 1, j)
            c = ring_idx(i + 1, j + 1)
            d = ring_idx(i, j + 1)
            lines.append("f %d//%d %d//%d %d//%d" % (a, a, b, b, d, d))
            lines.append("f %d//%d %d//%d %d//%d" % (b, b, c, c, d, d))

    # South pole cap: triangles connecting the last ring to the south pole.
    # Winding: ring_last[j+1], ring_last[j], south_pole  (CCW from outside)
    last_ring = num_lat - 1
    for j in range(num_lon):
        v0 = ring_idx(last_ring, j + 1)
        v1 = ring_idx(last_ring, j)
        v2 = south_pole
        lines.append("f %d//%d %d//%d %d//%d" % (v0, v0, v1, v1, v2, v2))

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.stdout.write(generate_sphere_obj(num_lon=16, num_lat=12, radius=1.0))
