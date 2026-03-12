-- Hollow pipe (tube) with all edges filleted.
-- Tests mixed convex and concave toroidal surfaces: filleting the circular edges
-- of a hollow pipe produces convex toroidal fillets on the outer edges and
-- concave toroidal fillets on the inner edges.
-- Outer diameter: 20 (r_out=10), inner diameter: 12 (r_in=6), height: 30.
-- Fillet radius: 1.5mm. Annulus width = 4mm, total fillet width = 3mm, no overlap.
-- Expected: 2 cylindrical surfaces (outer convex, inner concave) +
--   2 annular planar faces + 4 toroidal fillets (2 convex outer, 2 concave inner)
--   = 8 faces.
--
-- Analytical values (exact values from STEP file):
--   Outer torus: major R=8.5 (cyl_r - fillet_r), minor r=1.5
--   Inner torus: major R=7.5 (hole_r + fillet_r), minor r=1.5

local outer = cylinder(20, 30)    -- outer diameter 20 (r=10)
local inner = cylinder(12, 32)    -- inner diameter 12 (r=6), taller for clean cut
local pipe = difference(outer, inner)
local result = fillet_all(pipe, 1.5)
emit(result)
