-- 10mm cube. One corner at the origin (box places corner at origin).
-- Tests basic planar hypothesis fitting: 6 planar faces.
-- Expected: 6 planar hypotheses.
--
-- Analytical values:
--   Volume = 1000 mm^3
--   Surface area = 600 mm^2

local b = box(10, 10, 10)
emit(b)
