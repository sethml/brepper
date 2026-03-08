-- Sphere on top of a cylinder stalk (mushroom shape).
-- Cylinder: diameter 10 (r=5), height 20, centered in x/y, z=0 to z=20.
-- Sphere: diameter 16 (r=8), centered at (0,0,20) — sits on top of stalk.
-- The sphere radius > cylinder radius, so it bulges out beyond the stalk.
-- Tests sphere-cylinder boundary disambiguation.
-- Expected: 1 convex cylindrical surface (r=5) + 1 convex spherical surface (r=8) + 1 planar face (bottom cap).
--
-- Analytical values (intersection at z = 20 - sqrt(39)):
--   Volume = pi * [25*(20-sqrt(39)) + (8+sqrt(39))^2*(16-sqrt(39))/3] ≈ 3153 mm^3
--   Surface area = pi * [25 + 10*(20-sqrt(39)) + 16*(8+sqrt(39))] ≈ 1226 mm^2

local stalk = cylinder(10, 20)   -- diameter 10, height 20, z=0 to z=20
local ball = translate(sphere(16), 0, 0, 20)  -- diameter 16, center at (0,0,20)
local result = union(stalk, ball)
emit(result)
