-- Nosecone: truncated cone base with a tangent spherical cap at the tip.
-- Tests cone-sphere tangency: the conical and spherical surfaces meet tangent
-- along a circle (shared surface normals at the junction).
-- Expected: 1 convex conical face + 1 convex spherical face + 1 planar base.
--
-- Geometry:
--   Cone: half-angle alpha = 30 degrees (from axis), base radius 15, top (tangent circle) radius 5.
--   Cone height = (15 - 5) / tan(30) = 10*sqrt(3) ≈ 17.321 mm.
--   Sphere: tangent to cone at the top circle (r=5).
--     Sphere radius = 5 / cos(30) = 10/sqrt(3) ≈ 5.774 mm.
--     Sphere center at z = h + 5*tan(30) = 10*sqrt(3) + 5/sqrt(3) ≈ 20.207 mm.
--     Sphere top at z ≈ 25.981 mm.
--
-- Analytical values:
--   Cone lateral area = pi*(15+5)*sqrt((10*sqrt(3))^2 + (15-5)^2) = pi*20*sqrt(400) = 400*pi ≈ 1256.64 mm^2
--   Base cap area = pi*225 ≈ 706.86 mm^2
--   Spherical cap area = 2*pi*R_sphere*h_cap where h_cap = R_sphere - R_sphere*sin(30) = R_sphere*(1-0.5)
--     = 2*pi*(10/sqrt(3))*(10/(2*sqrt(3))) = 2*pi*100/(2*3) = 100*pi/3 ≈ 104.72 mm^2
--   Total surface area ≈ 2068.22 mm^2
--   Cone volume = pi/3 * h * (r1^2 + r1*r2 + r2^2) = pi/3 * 10*sqrt(3) * (225+75+25) = pi/3*10*sqrt(3)*325
--     = 3250*sqrt(3)*pi/3 ≈ 5890.49 mm^3
--   Spherical cap volume = pi*h_cap^2*(3*R-h_cap)/3 where h_cap = R*(1-sin(30)) = R/2
--     = pi*(R/2)^2*(3R-R/2)/3 = pi*R^2/4 * 5R/6 = 5*pi*R^3/24 = 5*pi*(10/sqrt(3))^3/24
--     = 5*pi*1000/(3*sqrt(3)*24) ≈ 100.72 mm^3
--   Total volume ≈ 5991.21 mm^3

local alpha = math.rad(30)               -- half-angle from axis
local r1 = 15                             -- base radius
local r2 = 5                              -- tangent circle radius
local h = (r1 - r2) / math.tan(alpha)    -- cone height = 10*sqrt(3)

-- Sphere tangent to cone at top circle
local R_sphere = r2 / math.cos(alpha)     -- sphere radius
local z_c = h + r2 * math.tan(alpha)      -- sphere center z

local c = cone(2*r1, 2*r2, h)
local s = translate(sphere(2*R_sphere), 0, 0, z_c)
local result = union(c, s)
emit(result)
