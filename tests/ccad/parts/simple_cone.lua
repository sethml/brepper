-- Truncated cone (frustum) centered in x/y, z from 0 to 30.
-- Tests basic conical surface fitting: 1 conical surface + 2 planar caps.
-- Expected: 1 convex conical hypothesis + 2 planar hypotheses.
-- Cone parameters: bottom d=20 (r=10), top d=10 (r=5), height=30.
-- Half-angle: atan((10-5)/30) = atan(1/6) ≈ 9.46°
--
-- Analytical values:
--   Volume = pi/3 * 30 * (100 + 50 + 25) = 1750*pi ≈ 5497.79 mm^3
--   Surface area = pi*15*sqrt(925) + pi*25 + pi*100 ≈ 1826.33 mm^2
--     (lateral: pi*15*sqrt(925) ≈ 1433.63, top cap: 25*pi ≈ 78.54, bottom cap: 100*pi ≈ 314.16)

local c = cone(20, 10, 30)  -- bottom diameter 20, top diameter 10, height 30
emit(c)
