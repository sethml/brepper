#!/usr/bin/env python3
"""Generate bad STL test files for brepper test suite.

Creates STL files in tests/bad/ that exercise specific failure modes:
- degenerate_face.stl: cube with one zero-area triangle
- non_manifold_edge.stl: three triangles sharing one edge
- inconsistent_winding.stl: cube with one triangle's winding reversed
- cube_shifted.stl: unit cube with one vertex shifted 10mm for --compare failure
"""

import os
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
BAD_DIR = os.path.join(PROJECT_DIR, "tests", "bad")


def write_ascii_stl(filename, triangles, solid_name="solid"):
    """Write an ASCII STL file.
    Each triangle is ((nx,ny,nz), (v1x,v1y,v1z), (v2x,v2y,v2z), (v3x,v3y,v3z)).
    """
    path = os.path.join(BAD_DIR, filename)
    with open(path, "w") as f:
        f.write(f"solid {solid_name}\n")
        for normal, v1, v2, v3 in triangles:
            f.write(f"  facet normal {normal[0]} {normal[1]} {normal[2]}\n")
            f.write("    outer loop\n")
            f.write(f"      vertex {v1[0]} {v1[1]} {v1[2]}\n")
            f.write(f"      vertex {v2[0]} {v2[1]} {v2[2]}\n")
            f.write(f"      vertex {v3[0]} {v3[1]} {v3[2]}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")
    print(f"  Created {path} ({len(triangles)} triangles)")


def cube_triangles(shift_vertex=None):
    """Generate 12 triangles for a unit cube [0,1]^3 with consistent outward normals.

    Winding follows the right-hand rule: vertices are ordered so that the
    cross product (v2-v1)x(v3-v1) points outward.

    Adjacent triangles share edges with opposite winding (manifold condition).

    If shift_vertex is (index, (dx,dy,dz)), vertex at that index is shifted.
    """
    V = [
        [0, 0, 0],  # 0
        [1, 0, 0],  # 1
        [1, 1, 0],  # 2
        [0, 1, 0],  # 3
        [0, 0, 1],  # 4
        [1, 0, 1],  # 5
        [1, 1, 1],  # 6
        [0, 1, 1],  # 7
    ]

    if shift_vertex is not None:
        idx, (dx, dy, dz) = shift_vertex
        V[idx] = [V[idx][0] + dx, V[idx][1] + dy, V[idx][2] + dz]

    # Convert to tuples
    V = [tuple(v) for v in V]

    # 12 triangles, 2 per face, outward normals.
    # For each face, the two triangles share a diagonal edge.
    # The winding ensures adj triangles traverse shared edges in opposite directions.
    triangles = [
        # Bottom face (z=0), outward normal = (0,0,-1)
        # Looking from below: vertices appear CCW as 0,3,2,1
        ((0, 0, -1), V[0], V[3], V[1]),  # T0: diagonal 0-3, edge 3-1, edge 1-0
        ((0, 0, -1), V[1], V[3], V[2]),  # T1: edge 1-3, diagonal 3-2 (wait...)
        # Top face (z=1), outward normal = (0,0,1)
        # Looking from above: vertices appear CCW as 4,5,6,7
        ((0, 0, 1), V[4], V[5], V[7]),   # T2
        ((0, 0, 1), V[5], V[6], V[7]),   # T3
        # Front face (y=0), outward normal = (0,-1,0)
        # Looking from front: x goes right, z goes up -> CCW = 0,1,5,4
        ((0, -1, 0), V[0], V[1], V[4]),  # T4
        ((0, -1, 0), V[1], V[5], V[4]),  # T5
        # Back face (y=1), outward normal = (0,1,0)
        # Looking from back: x goes left, z goes up -> CCW = 2,3,7,6
        ((0, 1, 0), V[2], V[3], V[6]),   # T6
        ((0, 1, 0), V[3], V[7], V[6]),   # T7
        # Left face (x=0), outward normal = (-1,0,0)
        # Looking from left: y goes right, z goes up -> CCW = 3,0,4,7
        ((-1, 0, 0), V[3], V[0], V[7]),  # T8
        ((-1, 0, 0), V[0], V[4], V[7]),  # T9
        # Right face (x=1), outward normal = (1,0,0)
        # Looking from right: y goes left, z goes up -> CCW = 1,2,6,5
        ((1, 0, 0), V[1], V[2], V[5]),   # T10
        ((1, 0, 0), V[2], V[6], V[5]),   # T11
    ]
    return triangles


def generate_degenerate_face():
    """Cube with one degenerate (zero-area) triangle.
    Replace the first triangle with one that has collinear vertices.
    """
    triangles = cube_triangles()
    # Replace T0 with a degenerate triangle: three collinear points on bottom face
    triangles[0] = ((0, 0, -1), (0, 0, 0), (0.5, 0, 0), (1, 0, 0))
    write_ascii_stl("degenerate_face.stl", triangles, "degenerate_cube")


def generate_non_manifold_edge():
    """Three triangles sharing one edge (a non-manifold configuration).
    Two triangles form a "tent" and a third shares their common edge.
    """
    triangles = [
        # Triangle 1: base edge (0,0,0)-(1,0,0), apex at (0.5, 1, 0)
        ((0, 0, 1), (0, 0, 0), (1, 0, 0), (0.5, 1, 0)),
        # Triangle 2: shares edge (0,0,0)-(1,0,0), apex at (0.5, -1, 0)
        ((0, 0, -1), (1, 0, 0), (0, 0, 0), (0.5, -1, 0)),
        # Triangle 3: ALSO shares edge (0,0,0)-(1,0,0), apex at (0.5, 0, 1)
        ((0, 0, 1), (0, 0, 0), (1, 0, 0), (0.5, 0, 1)),
    ]
    write_ascii_stl("non_manifold_edge.stl", triangles, "non_manifold")


def generate_inconsistent_winding():
    """Cube with one triangle having reversed winding (inconsistent orientation).
    This flips one triangle so its normal points inward instead of outward.
    """
    triangles = cube_triangles()
    # Reverse the winding of triangle T4 (front face, first triangle)
    # Original: ((0,-1,0), V[0], V[1], V[4])
    # Reversed: ((0,1,0), V[4], V[1], V[0])  -- swap v1 and v3
    n, v1, v2, v3 = triangles[4]
    triangles[4] = ((-n[0], -n[1], -n[2]), v3, v2, v1)
    write_ascii_stl("inconsistent_winding.stl", triangles, "inconsistent_cube")


def generate_cube_shifted():
    """Unit cube with vertex 6 (1,1,1) shifted to (5,5,5) for --compare failure.
    The shift must move the vertex away from ALL infinite planes of the cube
    (x=0, x=1, y=0, y=1, z=0, z=1), since ProjectPointOnSurf projects onto
    unbounded surfaces. The point (5,5,5) is at least 4mm from every plane.
    """
    triangles = cube_triangles(shift_vertex=(6, (4, 4, 4)))
    write_ascii_stl("cube_shifted.stl", triangles, "shifted_cube")

    # Copy the normal cube.step as reference for --compare
    src = os.path.join(PROJECT_DIR, "tests", "manual", "cube.step")
    dst = os.path.join(BAD_DIR, "cube_shifted.step")
    shutil.copy2(src, dst)
    print(f"  Copied {src} -> {dst}")


def main():
    os.makedirs(BAD_DIR, exist_ok=True)
    print("Generating bad test STL files...")
    generate_degenerate_face()
    generate_non_manifold_edge()
    generate_inconsistent_winding()
    generate_cube_shifted()
    print("Done.")


if __name__ == "__main__":
    main()
