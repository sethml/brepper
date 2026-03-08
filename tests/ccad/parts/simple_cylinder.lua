-- Simple cylinder with flat caps, centered in x/y, z from 0 to 30.
-- Tests basic cylindrical hypothesis fitting: 1 cylindrical surface + 2 planar caps.
-- Expected: 1 convex cylindrical hypothesis (radius 10) + 2 planar hypotheses.
--
-- Analytical values:
--   Volume = 3000*pi ≈ 9424.78 mm^3
--   Surface area = 800*pi ≈ 2513.27 mm^2

local c = cylinder(20, 30)  -- diameter 20, height 30
emit(c)
