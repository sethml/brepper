-- Block with a quarter-cylindrical notch cut from one corner.
-- The block has a corner at the origin and the cylinder is centered at the origin,
-- so only a quarter of the cylinder intersects the block corner.
-- Tests concave cylindrical surface with partial (quarter) arc + trimmed planar faces.
-- Expected: 6 trimmed planar faces + 1 concave quarter-cylindrical face.
--
-- Analytical values:
--   Volume = 18000 - 180*pi ≈ 17434.60 mm^3
--   Surface area = 3960 + 42*pi ≈ 4091.95 mm^2
--
-- Note: box(w,d,h) places a corner at the origin. cylinder(d,h) is centered in x/y.
-- So the cylinder at origin only overlaps the corner of the box.
local block = box(30, 30, 20)
local hole = cylinder(12, 22)  -- slightly taller to ensure clean through-cut
local result = difference(block, hole)
emit(result)
