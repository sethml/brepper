-- Stepped cylinder: two coaxial cylinders with different radii.
-- Bottom: diameter 20 (r=10), height 15. Top: diameter 14 (r=7), height 15 at z=15.
-- Tests that separate cylindrical hypotheses form for same-axis, different-radius surfaces.
-- Expected: 2 cylindrical surfaces + 3 planar surfaces (bottom cap, step ring, top cap).
--
-- Analytical values:
--   Volume = 2235*pi ≈ 7021.46 mm^3
--   Surface area = 710*pi ≈ 2230.53 mm^2

local bottom = cylinder(20, 15)
local top_part = translate(cylinder(14, 15), 0, 0, 15)
local result = union(bottom, top_part)
emit(result)
