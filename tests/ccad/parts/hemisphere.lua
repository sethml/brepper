-- Hemisphere: upper half of a sphere (z >= 0), flat base at z=0.
-- Tests partial convex spherical surface + plane-sphere boundary.
-- Expected: 1 convex spherical hypothesis (radius 10) + 1 planar hypothesis (flat bottom).
--
-- Analytical values:
--   Volume = 2000*pi/3 ≈ 2094.40 mm^3
--   Surface area = 300*pi ≈ 942.48 mm^2  (curved: 200*pi, flat: 100*pi)

local s = sphere(20)  -- diameter 20, radius 10, centered at origin
-- Cutting box covers z >= 0 region, extends well beyond the sphere
local half = translate(box(30, 30, 10), -15, -15, 0)
local result = intersection(s, half)
emit(result)
