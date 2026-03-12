-- Cone joined to cylinder: truncated cone base topped with a cylinder.
-- The cone's top diameter matches the cylinder's diameter for a smooth tangent junction.
-- Tests cone-cylinder tangency and mixed surface fitting.
-- Expected: 1 convex conical hypothesis + 1 convex cylindrical hypothesis + 2 planar caps.
-- Cone: bottom d=20 (r=10), top d=10 (r=5), height=15.
-- Cylinder: d=10 (r=5), height=15.
-- Half-angle: atan((10-5)/15) = atan(1/3) ≈ 18.43°
--
-- Analytical values:
--   Volume = pi/3 * 15 * (100+50+25) + pi * 25 * 15 = (875+375)*pi = 1250*pi ≈ 3926.99 mm^3
--   Surface area = pi*100 + pi*15*sqrt(250) + 2*pi*5*15 + pi*25 ≈ 1608.97 mm^2
--     (bottom cap: 100*pi ≈ 314.16, cone lateral: pi*15*sqrt(250) ≈ 745.04,
--      cylinder lateral: 150*pi ≈ 471.24, top cap: 25*pi ≈ 78.54)

local c = cone(20, 10, 15)
local cyl = translate(cylinder(10, 15), 0, 0, 15)
local result = union(c, cyl)
emit(result)
