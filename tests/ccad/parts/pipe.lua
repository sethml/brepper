-- Hollow tube (pipe section), centered in x/y, z from 0 to 30.
-- Tests convexity check: inner surface is concave, outer surface is convex.
-- Also tests annular planar end faces.
-- Expected: 2 cylindrical surfaces (convex outer r=10, concave inner r=7) + 2 annular planar faces.
--
-- Analytical values:
--   Volume = 1530*pi ≈ 4806.64 mm^3
--   Surface area = 1122*pi ≈ 3524.87 mm^2  (outer: 600*pi, inner: 420*pi, annuli: 102*pi)

local outer = cylinder(20, 30)   -- outer diameter 20 (r=10)
local inner = cylinder(14, 32)   -- inner diameter 14 (r=7), slightly taller for clean cut
local result = difference(outer, inner)
emit(result)
