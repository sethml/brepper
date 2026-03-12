-- Rectangular block with a centered conical through-hole (tapered bore).
-- Tests concave conical hypothesis + plane-cone boundary faces.
-- Expected: 6 planar faces (incl. 2 with different-diameter circular holes) + 1 concave conical face.
-- Cone parameters: bottom d=12 (r=6), top d=8 (r=4), block height=20.
-- Half-angle: atan((6-4)/20) = atan(0.1) ≈ 5.71°
--
-- Analytical values:
--   Volume = 18000 - pi/3 * 20 * (36 + 24 + 16) = 18000 - 1520*pi/3 ≈ 16408.90 mm^3
--   Surface area = 4200 - pi*52 + pi*10*sqrt(404) ≈ 4668.05 mm^2
--     (block outer: 4200, minus 2 holes: pi*(36+16)=52*pi ≈ 163.36,
--      conical bore lateral: pi*10*sqrt(404) ≈ 631.41)

local block = center_xyz(box(30, 30, 20))  -- centered at origin
-- cone(d1, d2, h): bottom d=12 at z=-11, top d=8 at z=+11
-- slightly taller than block to ensure clean through-cut
local hole = translate(cone(12, 8, 22), 0, 0, -11)
local result = difference(block, hole)
emit(result)
