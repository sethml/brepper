-- T-shape: two boxes forming a T when viewed from above (XY plane).
-- Base: 20x5x10 along x-axis. Stem: 5x10x10 extending in y from base center.
-- Both boxes have corners at the origin (box places corner at origin).
-- Tests planar fitting with coplanar faces and complex face adjacency.
--
-- Analytical values:
--   Volume = 1500 mm^3
--   Surface area = 1000 mm^2

local base = box(20, 5, 10)
local top_part = translate(box(5, 10, 10), 7.5, 5, 0)
local t = union(base, top_part)
emit(t)
