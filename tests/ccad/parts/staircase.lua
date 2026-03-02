-- Staircase: stepped block with 4 levels.
-- Tests planar fitting with many parallel planes at different heights.

local s1 = box(20, 10, 5)
local s2 = translate(box(15, 10, 5), 0, 0, 5)
local s3 = translate(box(10, 10, 5), 0, 0, 10)
local s4 = translate(box(5, 10, 5), 0, 0, 15)
local staircase = union(union(union(s1, s2), s3), s4)
emit(staircase)
