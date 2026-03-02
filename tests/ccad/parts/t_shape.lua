-- T-shape: two boxes forming a T.
-- Tests planar fitting with coplanar faces and complex face adjacency.

local base = box(20, 5, 10)
local top_part = translate(box(5, 10, 10), 7.5, 5, 0)
local t = union(base, top_part)
emit(t)
