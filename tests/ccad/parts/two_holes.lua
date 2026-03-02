-- Block with two cylindrical through-holes of different diameters.
-- Tests multiple independent concave cylindrical hypotheses at different locations.
-- Expected: 6 planar faces + 2 concave cylindrical faces (different radii).

local block = box(40, 20, 20)
local hole1 = translate(cylinder(10, 22), -10, 0, 0)
local hole2 = translate(cylinder(6, 22), 10, 0, 0)
local result = difference(difference(block, hole1), hole2)
emit(result)
