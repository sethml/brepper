-- Stepped cylinder: two coaxial cylinders with different radii.
-- Tests that separate cylindrical hypotheses form for same-axis, different-radius surfaces.
-- Expected: 2 cylindrical surfaces + 3 planar surfaces (bottom cap, step ring, top cap).

local bottom = cylinder(20, 15)
local top_part = translate(cylinder(14, 15), 0, 0, 15)
local result = union(bottom, top_part)
emit(result)
