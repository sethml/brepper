-- Rectangular block with a centered cylindrical through-hole.
-- Tests concave cylindrical hypothesis + plane-cylinder boundary faces.
-- Expected: 6 planar faces (incl. 2 annular) + 1 concave cylindrical face.
--
-- Analytical values:
--   Volume = 18000 - 720*pi ≈ 15737.17 mm^3
--   Surface area = 4200 + 168*pi ≈ 4727.79 mm^2
--
-- box(w,d,h) places a corner at the origin, so we center it.
-- cylinder(d,h) is already centered in x/y.
local block = center_xyz(box(30, 30, 20))  -- centered at origin
local hole = translate(cylinder(12, 22), 0, 0, -11)  -- slightly taller to ensure clean through-cut
local result = difference(block, hole)
emit(result)
