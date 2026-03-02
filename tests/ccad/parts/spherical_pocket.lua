-- Block with hemispherical pocket carved from the top.
-- Tests concave spherical surface detection.
-- Expected: 6 planar hypotheses (including annular top face) + 1 concave spherical hypothesis (radius 8).

local block = box(30, 30, 20)
-- Sphere centered at top-center of block, radius 8
local pocket = translate(sphere(16), 15, 15, 20)
local result = difference(block, pocket)
emit(result)
