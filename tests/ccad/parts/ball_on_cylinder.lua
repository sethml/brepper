-- Sphere on top of a cylinder stalk.
-- Tests sphere-cylinder boundary disambiguation.
-- Expected: 1 convex cylindrical hypothesis (radius 5) + 1 convex spherical hypothesis (radius 8) + 1 planar hypothesis (bottom cap).

local stalk = cylinder(10, 20)   -- diameter 10, height 20, z=0 to z=20
local ball = translate(sphere(16), 0, 0, 20)  -- diameter 16, center at (0,0,20)
local result = union(stalk, ball)
emit(result)
