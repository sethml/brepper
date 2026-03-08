-- Simple sphere, centered at origin.
-- Tests basic convex spherical hypothesis fitting: 1 spherical surface.
-- Expected: 1 convex spherical hypothesis with radius 10.
--
-- Analytical values:
--   Volume = 4000*pi/3 ≈ 4188.79 mm^3
--   Surface area = 400*pi ≈ 1256.64 mm^2

local s = sphere(20)  -- diameter 20, radius 10, centered at origin
emit(s)
