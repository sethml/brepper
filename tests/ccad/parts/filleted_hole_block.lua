-- Block with cylindrical through-hole, hole edges filleted.
-- Tests concave toroidal surface detection: filleting the circular edges
-- where a cylindrical bore meets the top/bottom planar faces creates
-- concave toroidal fillet surfaces inside the bore.
-- Block: 30×30×20, hole diameter: 14 (r=7), fillet radius: 2mm.
-- Expected: 6 planar faces (top/bottom trimmed by hole) + 1 concave cylindrical
--   surface (hole wall) + 2 concave toroidal fillets = 9 faces.
--
-- Analytical values (exact values from STEP file):
--   Torus parameters: major R=9 (hole_r + fillet_r), minor r=2 (fillet_r)

local block = box(30, 30, 20)
local hole = translate(cylinder(14, 22), 15, 15, -1)  -- taller for clean cut
local part = difference(block, hole)
-- Select only the circular edges (hole-plane intersections)
local circ_edges = edges(part):geom("circle"):collect()
local result = fillet(part, circ_edges, 2)
emit(result)
