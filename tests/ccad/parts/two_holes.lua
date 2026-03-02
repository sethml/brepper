-- Block with two cylindrical through-holes of different diameters.
-- Tests multiple independent concave cylindrical hypotheses at different locations.
-- Expected: 6 planar faces + 2 concave cylindrical faces (different radii).
-- box() places origin at corner, so we center the block first.

local block = translate(box(40, 20, 20), -20, -10, -10)
local hole1 = translate(cylinder(10, 22), -10, 0, 0)
local hole2 = translate(cylinder(6, 22), 10, 0, 0)
local result = difference(difference(block, hole1), hole2)
emit(result)
