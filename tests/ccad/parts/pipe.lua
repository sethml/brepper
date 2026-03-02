-- Hollow tube (pipe section).
-- Tests convexity check: inner surface is concave, outer surface is convex.
-- Also tests annular planar end faces.
-- Expected: 2 cylindrical surfaces (inner + outer) + 2 annular planar faces.

local outer = cylinder(20, 30)
local inner = cylinder(14, 32)  -- slightly taller for clean cut
local result = difference(outer, inner)
emit(result)
