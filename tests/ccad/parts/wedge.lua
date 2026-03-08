-- Wedge with non-axis-aligned planar faces.
-- wedge(dx, dy, dz, ltx) creates a box that narrows: bottom is dx*dy, top is ltx*dy.
-- The top face is centered in x (offset by (dx-ltx)/2 = 2.5mm from each side).
-- Corner at origin (wedge places corner at origin).
-- Tests planar fitting with angled surfaces.
--
-- Analytical values (from OCCT BRepGProp):
--   Volume = 750 mm^3
--   Surface area ≈ 511.80 mm^2

local w = wedge(10, 10, 10, 5)
emit(w)
