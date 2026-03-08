-- Staircase: 4 descending steps from left (x=0) to right (x=20).
-- Each step is 5mm wide (in x), 10mm deep (y), and 5mm tall (z).
-- Tallest step is at x=0 (height 20mm), shortest at x=15-20 (height 5mm).
-- All boxes have corners at the origin (box places corner at origin).
-- Tests planar fitting with many parallel planes at different heights.
--
-- Analytical values:
--   Volume = 2500 mm^3
--   Surface area = 1300 mm^2

local s1 = box(20, 10, 5)
local s2 = translate(box(15, 10, 5), 0, 0, 5)
local s3 = translate(box(10, 10, 5), 0, 0, 10)
local s4 = translate(box(5, 10, 5), 0, 0, 15)
local staircase = union(union(union(s1, s2), s3), s4)
emit(staircase)
