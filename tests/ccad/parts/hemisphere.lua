-- Hemisphere: top half of a sphere on a flat base.
-- Tests partial convex spherical surface + plane-sphere boundary.
-- Expected: 1 convex spherical hypothesis (radius 10) + 1 planar hypothesis (flat bottom).

local s = sphere(20)  -- diameter 20, radius 10, centered at origin
-- Cut box covering top half only (z >= 0)
local half = translate(box(30, 30, 10), -15, -15, 0)
local result = intersection(s, half)
emit(result)
