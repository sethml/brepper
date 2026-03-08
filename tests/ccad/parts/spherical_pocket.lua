-- Block with hemispherical pocket carved from the top.
-- The block has a corner at the origin. The sphere is centered at the
-- top-center of the block (15, 15, 20), creating a concave hemispherical
-- pocket in the top face.
-- Tests concave spherical surface detection.
-- Expected: 6 planar faces (including annular top face) + 1 concave spherical face (radius 8).
--
-- Analytical values:
--   Volume = 18000 - 1024*pi/3 ≈ 16928.05 mm^3
--   Surface area = 4200 + 64*pi ≈ 4401.06 mm^2

local block = box(30, 30, 20)
-- Sphere centered at top-center of block, radius 8
local pocket = translate(sphere(16), 15, 15, 20)
local result = difference(block, pocket)
emit(result)
