-- Simple sphere.
-- Tests basic convex spherical hypothesis fitting: 1 spherical surface.
-- Expected: 1 convex spherical hypothesis with radius 10.

local s = sphere(20)  -- diameter 20, radius 10, centered at origin
emit(s)
