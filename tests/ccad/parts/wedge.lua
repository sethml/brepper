-- Wedge (triangular prism) with non-axis-aligned planar faces.
-- wedge(dx, dy, dz, ltx) creates a box that narrows from dx to ltx along y.
-- With ltx = dx/2, the top edge is a single line, forming a triangular prism.
-- Corner at origin (wedge places corner at origin).
-- Tests planar fitting with angled surfaces.
--
-- Analytical values:
--   Volume = 500 mm^3
--   Surface area = 200 + 100*sqrt(5) ≈ 423.61 mm^2

local w = wedge(10, 10, 10, 5)
emit(w)
