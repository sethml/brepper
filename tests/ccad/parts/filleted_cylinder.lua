-- Cylinder with all edges filleted.
-- Tests convex toroidal surface detection: filleting the circular edges of a
-- cylinder creates toroidal fillet surfaces where the cylindrical wall meets
-- the planar end caps. This is the simplest model with toroidal surfaces.
-- Cylinder: d=20, h=30, fillet radius: 3mm.
-- Expected: 1 cylindrical surface + 2 planar end caps + 2 convex toroidal fillets = 5 faces.
--
-- Analytical values:
--   Cylindrical surface: 2*pi*10*24 ≈ 1507.96 mm^2
--   2 end caps: 2*pi*7^2 ≈ 307.88 mm^2
--   2 quarter-torus fillets (R=7, r=3): 2*pi^2*7*3 ≈ 414.52 mm^2
--   Total surface area ≈ 2230.36 mm^2
--   Volume (Pappus) ≈ 9198.33 mm^3

local c = cylinder(20, 30)   -- diameter 20, height 30
local result = fillet_all(c, 3)
emit(result)
