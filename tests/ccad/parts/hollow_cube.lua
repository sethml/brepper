-- 20mm cube with a fully-enclosed 10mm cubic void at its center.
-- Tests multi-shell solid construction: the mesh has 2 closed shells
-- (outer cube surface + inner void surface), and the STEP output must be
-- ONE solid with an internal cavity, not two nested solids.
-- Expected: 12 planar surfaces (6 outer + 6 inner void).
--
-- Analytical values:
--   Volume = 8000 - 1000 = 7000 mm^3
--   Surface area = 2400 + 600 = 3000 mm^2

local outer = box(20, 20, 20)
local void = translate(box(10, 10, 10), 5, 5, 5)
local result = difference(outer, void)
emit(result)
