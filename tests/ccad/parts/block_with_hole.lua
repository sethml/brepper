-- Rectangular block with a cylindrical through-hole.
-- Tests concave cylindrical hypothesis + plane-cylinder boundary faces.
-- Expected: 6 planar faces + 1 concave cylindrical face.

local block = box(30, 30, 20)
local hole = cylinder(12, 22)  -- slightly taller to ensure clean through-cut
local result = difference(block, hole)
emit(result)
