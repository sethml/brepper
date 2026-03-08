-- Block with two cylindrical through-holes of different diameters.
-- The block is centered at the origin. Holes are at x=-10 and x=+10.
-- Tests multiple independent concave cylindrical hypotheses at different locations.
-- Expected: 6 planar faces + 2 concave cylindrical faces (r=5 and r=3).
--
-- Analytical values:
--   Volume = 16000 - 680*pi ≈ 13863.72 mm^3
--   Surface area = 4000 + 252*pi ≈ 4791.68 mm^2
--
-- box() places corner at origin; cylinder() is centered in x/y. We translate
-- both to center the block and position holes symmetrically.

local block = translate(box(40, 20, 20), -20, -10, -10)
local hole1 = translate(cylinder(10, 22), -10, 0, -11)
local hole2 = translate(cylinder(6, 22), 10, 0, -11)
local result = difference(difference(block, hole1), hole2)
emit(result)
