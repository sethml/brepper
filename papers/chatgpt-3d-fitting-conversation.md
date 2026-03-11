You: 

Find me an algorithm to fit a 3d cylinder to a set of 3d points minimizing sum of error squared.

ChatGPT: 

The standard formulation is a **nonlinear least-squares fit of a cylinder surface**. A cylinder can be parameterized by:

-   axis direction **a** (unit vector)
    
-   a point on the axis **q**
    
-   radius **r**
    

For each point $p_i$, the distance to the cylinder surface is:

$$
d_i = \left\| (p_i - q) \times a \right\| - r
$$

because $\|(p_i-q)\times a\|$ is the perpendicular distance from the point to the cylinder axis. [rdrr.io](https://rdrr.io/github/tiagodc/TreeLS/man/sgt.irls.cylinder.html?utm_source=chatgpt.com)

You then minimize

$$
E = \sum_i d_i^2
$$

This is a **6-parameter nonlinear least-squares problem**:

-   2 DOF for axis direction (unit vector constraint)
    
-   3 DOF for axis point
    
-   1 DOF for radius
    

---

# Practical algorithm (works well in practice)

## 1\. Estimate cylinder axis (good initialization)

Compute the covariance of the points:

$$
C = \frac{1}{N}\sum (p_i-\bar p)(p_i-\bar p)^T
$$

Let $v_1,v_2,v_3$ be eigenvectors.

For cylindrical point clouds:

-   **axis ≈ eigenvector with largest eigenvalue**
    

because the data spreads most along the axis.

```
center = mean(points)  
C = covariance(points)  
axis = eigenvector\_largest(C)  
axis = normalize(axis)
```

---

## 2\. Project points to plane perpendicular to axis

Choose an orthonormal basis $u,v$ spanning the plane perpendicular to **axis**.

For each point:

$$
x_i = (p_i - center)\cdot u
$$
 
$$
y_i = (p_i - center)\cdot v
$$

Now the problem becomes **2D circle fitting**.

---

## 3\. Fit circle in that plane

Use a linear least-squares circle fit (Taubin or Pratt).

Solve:

$$
x^2 + y^2 + Ax + By + C = 0
$$

Center:

$$
(x_c,y_c)=(-A/2,-B/2)
$$

Radius:

$$
r = \sqrt{x_c^2 + y_c^2 - C}
$$

Axis point:

$$
q = center + x_c u + y_c v
$$

---

## 4\. Refine with nonlinear least squares

Use **Levenberg–Marquardt** or **Nelder–Mead** to minimize

$$
E(a,q,r)=\sum_i \left(\|(p_i-q)\times a\|-r\right)^2
$$

The optimization variables:

```
axis (2 parameters)  
axis\_point (3 parameters)  
radius (1 parameter)
```

Libraries like **Ceres**, **Eigen LM**, or **scipy.optimize** work well.

---

# Pseudocode

```
C++

axis = PCA\_axis(points)  
center = mean(points)  
  
(u,v) = orthonormal\_basis(axis)  
  
for p in points:  
    x = dot(p\-center, u)  
    y = dot(p\-center, v)  
    circle\_points.append((x,y))  
  
(cx,cy,r) = least\_squares\_circle(circle\_points)  
  
q = center + cx\*u + cy\*v  
  
optimize axis,q,r to minimize:  
    sum((norm(cross(p\-q,axis))\-r)^2)
```

---

# Robust variant (for real scans)

If your points contain outliers:

1.  **RANSAC cylinder fit**
    
2.  Then **least-squares refinement**
    

This is common in point-cloud processing.

---

# If you're implementing this in C++ (given your OCCT work)

The most stable pipeline is:

1.  PCA axis estimate
    
2.  Pratt circle fit
    
3.  Levenberg–Marquardt refinement
    

This converges extremely reliably.

---

💡 If you'd like, I can also show you a **very compact closed-form method (Eberly / David Eberly algorithm)** that avoids iterative optimization for most of the solve. It's one of the nicest cylinder fits and popular in geometry kernels.

You: 

Is there a rust library that implements this? Or OCCT function?

ChatGPT: 

Short answer: **no widely-used Rust crate directly implements cylinder fitting yet**, and **OCCT also does not provide a built-in “fit cylinder to points” function**. You generally have to implement it yourself or use a point-cloud library.

Below is the practical landscape.

---

# Rust libraries

Right now the Rust ecosystem for **geometric primitive fitting** is still pretty thin.

### 1\. General point-cloud crates (no cylinder fitting built in)

Examples:

-   `threecrate` — point cloud + mesh processing framework [lib.rs](https://lib.rs/crates/threecrate?utm_source=chatgpt.com)
    
-   `epoint` — point cloud storage/processing crate [docs.rs](https://docs.rs/epoint/latest/epoint/?utm_source=chatgpt.com)
    

These give you data structures but **not primitive fitting algorithms**.

Typical Rust approach today:

```
nalgebra  
\+ levenberg-marquardt  
\+ your own cylinder residual
```

or

```
argmin  
levenberg-marquardt
```

for nonlinear least squares.

If you're already using `nalgebra`, this ends up being ~150–200 LOC.

---

### 2\. Rust wrappers around C++ point-cloud libraries

You could bind to libraries that already do primitive fitting.

The most common one is **Point Cloud Library (PCL)**, which contains many point-cloud algorithms including RANSAC-based primitive fitting and segmentation. [en.wikipedia.org](https://en.wikipedia.org/wiki/Point_Cloud_Library?utm_source=chatgpt.com)

PCL has:

```
pcl::SACSegmentationFromNormals  
pcl::SampleConsensusModelCylinder
```

but you’d need a Rust wrapper (none widely maintained right now).

---

# OCCT (OpenCascade)

OCCT is great at **analytic geometry once parameters are known**, but it **does not include primitive fitting algorithms** for point clouds.

You can create cylinders once you know the axis + radius:

```
gp\_Ax2  
gp\_Cylinder  
BRepPrimAPI\_MakeCylinder
```

because OCCT provides constructors for cylinders from axis and radius. [dev.opencascade.org](https://dev.opencascade.org/doc/occt-6.7.0/overview/html/user_guides__modeling_data.html?utm_source=chatgpt.com)

But there is **no API equivalent to**:

```
GeomAPI\_PointsToCylinder
```

Only things like:

```
GeomAPI\_PointsToBSpline  
GeomAPI\_PointsToBSplineSurface
```

for spline fitting.

So the typical OCCT workflow is:

```
points  
 -> external fitting algorithm  
 -> gp\_Ax2 axis  
 -> gp\_Cylinder  
 -> build surface/solid
```

---

# What people usually do in CAD kernels

Typical pipeline:

```
1\. PCA -> axis guess  
2\. 2D circle fit in orthogonal plane  
3\. LM refinement  
4\. build gp\_Cylinder
```

Then use it inside OCCT.

FreeCAD and reverse-engineering tools do something similar internally.

---

# If you want a Rust-native approach

The easiest stack is:

```
nalgebra  
levenberg-marquardt  
kdtree (optional)
```

Define residual:

```
ri = |(pi - q) × a| - r
```

Minimize:

```
sum(ri²)
```

---

✅ **My recommendation for your case (given your OCCT work):**

Implement the fit in Rust/C++ yourself using:

```
PCA (axis)  
Pratt circle fit  
LM refine
```

Then construct:

```
gp\_Ax2 axis  
gp\_Cylinder
```

It's surprisingly small code.

---

💡 If you want, I can also give you a **very robust industrial cylinder-fit implementation (~80 lines of C++)** that I’ve seen used in CAD kernels and reverse-engineering tools. It avoids most of the LM pitfalls and converges fast.

You: 

Summarize the rust-native approach in a form suitable for ingestion by an LLM coding agent.

ChatGPT: 

Below is a **compact spec suitable for an LLM coding agent**. It assumes Rust and common linear-algebra crates.

---

# Task

Fit a **3D cylinder** to a set of 3D points by minimizing **sum of squared radial residuals**.

Cylinder parameters:

```
axis\_point q: Vec3  
axis\_dir   a: Vec3 (unit vector)  
radius     r: f64
```

Residual for point `p`:

```
ri = || (p - q) × a || - r
```

Objective:

```
minimize Σ ri²
```

---

# Dependencies

Use common Rust math libraries:

```
nalgebra  
levenberg-marquardt (or argmin)
```

Required types:

```
Vector3<f64>  
Matrix3<f64>
```

---

# Algorithm Overview

```
1\. Compute PCA of points → initial cylinder axis  
2\. Project points into plane orthogonal to axis  
3\. Fit 2D circle to projected points  
4\. Convert circle center back to 3D axis point  
5\. Run Levenberg–Marquardt refinement on (axis, point, radius)
```

---

# Step 1 — PCA axis estimate

Compute centroid:

```
c = mean(points)
```

Compute covariance matrix:

```
C = Σ (p - c)(p - c)^T / N
```

Eigen-decompose `C`.

Initial axis direction:

```
a = eigenvector with largest eigenvalue  
a = normalize(a)
```

---

# Step 2 — Build orthonormal plane basis

Create two vectors orthogonal to `a`.

Example:

```
u = normalize(any\_vector\_not\_parallel\_to(a) × a)  
v = a × u
```

Now `{u, v, a}` is an orthonormal frame.

---

# Step 3 — Project points into plane

For each point:

```
d = p - c  
x = dot(d, u)  
y = dot(d, v)
```

Collect `(x, y)`.

---

# Step 4 — 2D circle least squares

Solve circle equation:

```
x² + y² + A x + B y + C = 0
```

Linear least squares.

Circle center:

```
cx = -A/2  
cy = -B/2
```

Radius:

```
r = sqrt(cx² + cy² - C)
```

Axis point in 3D:

```
q = c + cx\*u + cy\*v
```

---

# Step 5 — Nonlinear refinement

Optimize parameters:

```
θ = \[axis\_direction(2 params), axis\_point(3), radius(1)\]
```

Represent axis direction using **spherical angles** or normalize after each update.

Residual function:

```
for p in points:  
    d = p - q  
    radial = norm(cross(d, a))  
    ri = radial - r
```

Return vector `[ri]`.

Minimize:

```
Σ ri²
```

Use **Levenberg–Marquardt**.

---

# Termination criteria

Typical:

```
max\_iterations = 50  
gradient\_tol   = 1e-10  
parameter\_tol  = 1e-10
```

---

# Output

Return:

```
struct CylinderFit {  
    axis\_point: Vector3<f64>,  
    axis\_dir:   Vector3<f64>, // normalized  
    radius:     f64,  
    rms\_error:  f64  
}
```

---

# Numerical notes

-   Require at least **6 points**
    
-   Normalize axis after updates
    
-   Reject degenerate cases:
    
    -   very small radius
        
    -   nearly planar point sets
        

---

# Complexity

```
PCA: O(N)  
Circle fit: O(N)  
LM iterations: O(kN)   (k ≈ 10–30)
```

---

If useful, I can also produce a **ready-to-paste Rust implementation (~120 lines)** using `nalgebra` + `levenberg-marquardt`. It works well for CAD / scan fitting.

You: 

Restate that in raw markdown.

ChatGPT: 

```
Markdown

\# Task  
  
Fit a \*\*3D cylinder\*\* to a set of 3D points by minimizing \*\*sum of squared radial residuals\*\*.  
  
Cylinder parameters:  
  
\- axis\_point \`q: Vec3\`  
\- axis\_dir \`a: Vec3\` (unit vector)  
\- radius \`r: f64\`  
  
Residual for point \`p\`:
```

ri = || (p - q) × a || - r

```
Objective:
```

minimize Σ ri²

```
\---  
  
\# Dependencies  
  
Use common Rust math libraries:  
  
\- \`nalgebra\`  
\- \`levenberg-marquardt\` (or \`argmin\`)  
  
Required types:  
  
\- \`Vector3<f64>\`  
\- \`Matrix3<f64>\`  
  
\---  
  
\# Algorithm Overview
```

1.  Compute PCA of points → initial cylinder axis
    
2.  Project points into plane orthogonal to axis
    
3.  Fit 2D circle to projected points
    
4.  Convert circle center back to 3D axis point
    
5.  Run Levenberg–Marquardt refinement on (axis, point, radius)
    

```
\---  
  
\# Step 1 — PCA axis estimate  
  
Compute centroid:
```

c = mean(points)

```
Compute covariance matrix:
```

C = Σ (p - c)(p - c)^T / N

```
Eigen-decompose \`C\`.  
  
Initial axis direction:
```

a = eigenvector with largest eigenvalue  
a = normalize(a)

```
\---  
  
\# Step 2 — Build orthonormal plane basis  
  
Create two vectors orthogonal to \`a\`.  
  
Example:
```

u = normalize(any\_vector\_not\_parallel\_to(a) × a)  
v = a × u

```
Now \`{u, v, a}\` is an orthonormal frame.  
  
\---  
  
\# Step 3 — Project points into plane  
  
For each point:
```

d = p - c  
x = dot(d, u)  
y = dot(d, v)

```
Collect \`(x, y)\`.  
  
\---  
  
\# Step 4 — 2D circle least squares  
  
Solve circle equation:
```

x² + y² + A x + B y + C = 0

```
Using linear least squares.  
  
Circle center:
```

cx = -A/2  
cy = -B/2

```
Radius:
```

r = sqrt(cx² + cy² - C)

```
Axis point in 3D:
```

q = c + cx*u + cy*v

```
\---  
  
\# Step 5 — Nonlinear refinement  
  
Optimize parameters:
```

θ = \[axis\_direction(2 params), axis\_point(3), radius(1)\]

```
Represent axis direction using spherical angles or normalize after each update.  
  
Residual function:
```

for p in points:  
d = p - q  
radial = norm(cross(d, a))  
ri = radial - r

```
Return vector \`\[ri\]\`.  
  
Minimize:
```

Σ ri²

```
Use \*\*Levenberg–Marquardt\*\*.  
  
\---  
  
\# Termination criteria  
  
Typical values:
```

max\_iterations = 50  
gradient\_tol = 1e-10  
parameter\_tol = 1e-10

```
\---  
  
\# Output
```

struct CylinderFit {  
axis\_point: Vector3<f64>,  
axis\_dir: Vector3<f64>, // normalized  
radius: f64,  
rms\_error: f64  
}

```
\---  
  
\# Numerical notes  
  
\- Require at least \*\*6 points\*\*  
\- Normalize axis after updates  
\- Reject degenerate cases:  
  - very small radius  
  - nearly planar point sets  
  
\---  
  
\# Complexity
```

PCA: O(N)  
Circle fit: O(N)  
LM iterations: O(kN) (k ≈ 10–30)

You: 

I want to move on to toroid fits, where I'm fitting a toroidal fragment to a triangulated mesh - for example, a filleted edge that goes around a rounded corner. What algorithm should i use? What's the minimum number of points required? Does it make sense to use the same algorithm to detect toroidal fits and cylindrical fits, with the cylinder as a degenerate infinite-radius torus? Does it make sense to use the same algorithm to detect toroidal fits and spherical fits, with the sphere as a degenerate zero-radius toroid?

ChatGPT: 

Toroid fitting is a bit trickier than cylinders because the geometry has **two radii and a curved centerline**, but there’s a well-known approach that works well for CAD-like meshes such as fillets.

I’ll break it into the practical pieces.

---

# 1\. Toroid parameterization

A torus is defined by:

-   **Center** `C`
    
-   **Axis direction** `a` (unit vector)
    
-   **Major radius** `R` (distance from center to tube centerline)
    
-   **Minor radius** `r` (tube radius)
    

Implicit equation for point `p`:

```
let d = p - C  
let axial = dot(d, a)  
let radial\_vec = d - axial \* a  
let radial = norm(radial\_vec)  
  
ri = sqrt((radial - R)^2 + axial^2) - r
```

Residual:

```
ri
```

Minimize:

```
Σ ri²
```

This is the **standard torus least-squares residual**.

Unknowns:

```
C (3)  
a (2 dof)  
R  
r
```

Total: **7 parameters**

---

# 2\. Minimum number of points

Each point contributes **1 equation**.

Parameters = 7  
So the theoretical minimum is:

```
7 points
```

But this is **very unstable**.

Practical minimums:

| Situation | Points |
| --- | --- |
| minimal solve | 7 |
| stable fit | ~20 |
| robust mesh fitting | 50–200 |

For a triangulated fillet patch you usually have plenty.

---

# 3\. Good algorithm (recommended)

Use a **two-stage approach**:

### Stage 1 — estimate torus axis

Compute PCA of points.

For most torus fragments:

-   smallest eigenvalue → **axis direction**
    

This works because the torus surface bends around the axis.

---

### Stage 2 — reduce to 2D circle fitting

Project each point to cylindrical coordinates around the axis.

For each point:

```
d = p - C\_guess  
axial = dot(d, a)  
  
radial = norm(d - axial\*a)
```

Then treat:

```
(radial, axial)
```

as 2D coordinates.

The torus cross-section is a **circle**:

```
(radial - R)^2 + axial^2 = r^2
```

So fit a **2D circle** to points:

```
x = radial  
y = axial
```

This gives:

```
center = (R, 0)  
radius = r
```

---

### Stage 3 — nonlinear refinement

Then run **Levenberg–Marquardt** optimizing:

```
C  
a  
R  
r
```

Residual:

```
sqrt((radial - R)^2 + axial^2) - r
```

This converges very reliably.

---

# 4\. Important practical improvement

For **fillet patches**, the torus center `C` usually lies near the intersection of two offset surfaces.

But you can initialize it simply as:

```
C = centroid(points)
```

The LM stage will correct it.

---

# 5\. Degenerate cases

Now your question about **unifying primitives** is interesting.

## Cylinder as infinite-radius torus?

Mathematically:

```
R → ∞
```

But numerically this is **terrible**.

Problems:

-   parameters blow up
    
-   optimization becomes ill-conditioned
    
-   gradients vanish
    

So **do NOT model cylinders as torus fits**.

Use a separate cylinder model.

---

## Sphere as zero-major-radius torus?

If:

```
R = 0
```

Then the torus equation becomes:

```
sqrt(radial² + axial²) = r
```

which is a **sphere**.

So mathematically:

```
sphere = torus with R = 0
```

But again numerically it's not ideal.

Reason:

-   torus axis becomes undefined
    
-   solver degeneracy
    

So in practice:

| primitive | separate fit |
| --- | --- |
| sphere | yes |
| cylinder | yes |
| torus | yes |

---

# 6\. Practical CAD primitive detection strategy

Most reverse-engineering tools do:

```
1\. plane  
2\. sphere  
3\. cylinder  
4\. cone  
5\. torus
```

Each with its own solver.

Then pick the primitive with **lowest residual error**.

This is far more stable than trying to unify them.

---

# 7\. Detecting toroidal fillets

For mesh fragments that look like fillets:

Heuristics:

```
principal curvature signs same  
curvature varies strongly in one direction  
surface normal rotates around an axis
```

That usually isolates torus candidates before fitting.

---

# 8\. Complexity

For N points:

```
PCA: O(N)  
circle fit: O(N)  
LM: O(kN) (k≈20)
```

Very cheap.

---

# 9\. Summary

**Recommended torus fitting pipeline**

```
1\. PCA → axis  
2\. Convert points to (radial, axial)  
3\. 2D circle fit → R,r  
4\. Initial center = centroid  
5\. LM optimize {C,a,R,r}
```

Minimum points:

```
7 theoretical  
20+ recommended
```

Primitive degeneracy:

| Primitive | Torus degenerate? | Use same solver? |
| --- | --- | --- |
| sphere | yes (R=0) | not recommended |
| cylinder | yes (R→∞) | not recommended |

---

💡 If you're doing **mesh-to-CAD reconstruction**, I can also show you a **much more powerful trick used in reverse-engineering tools**: detecting toroidal fillets by intersecting **normal lines** instead of fitting the surface directly.

That method is *dramatically more robust for small torus patches* (like rounded corners).

You: 

Yes, explain that trick, with enough detail to implement it.

ChatGPT: 

This trick is widely used in **reverse-engineering pipelines** because toroidal fillets are often **small fragments of a torus**, and direct nonlinear fitting can be unstable. The key observation is about **surface normals**.

For a torus, **every surface normal intersects the torus tube centerline** (the circle of radius $R$). That gives a way to estimate the torus geometry from **normal intersections**, which is much more stable than fitting the implicit torus equation directly.

Below is a practical algorithm you can implement.

---

# Key geometric fact

For a torus with:

-   center $C$
    
-   axis $a$
    
-   major radius $R$
    
-   minor radius $r$
    

Each point $p$ on the surface has a normal $n$ pointing toward the **center of the tube cross-section circle**.

That tube center lies at:

$$
k = p + r\,n
$$

All such points $k$ lie on the **major circle**:

-   plane orthogonal to axis $a$
    
-   center $C$
    
-   radius $R$
    

So if we estimate those $k$ points, we can **fit a circle to them**.

The difficulty: **we don’t know $r$** yet.

But we can solve for it.

---

# Algorithm overview

```
1\. Collect mesh vertices and normals  
2\. Estimate torus axis with PCA on normals  
3\. Estimate tube radius r using normal intersections  
4\. Compute tube-center points k = p + r\*n  
5\. Fit a 3D circle to the k points → major circle  
6\. Derive torus parameters  
7\. Optionally refine with nonlinear least squares
```

---

# Step 1 — Input

For each mesh vertex:

```
p\_i : Vector3  
n\_i : Vector3 (unit normal)
```

Normals should be oriented consistently.

---

# Step 2 — Estimate torus axis

Normals of a torus **rotate around the axis**.

Compute PCA of normals:

```
C = Σ n\_i n\_i^T
```

Eigenvectors:

```
smallest eigenvalue → axis direction a
```

Normalize `a`.

---

# Step 3 — Estimate minor radius r

Two nearby points on the torus have normals that intersect near the tube center.

For two points:

```
L1: p1 + t n1  
L2: p2 + s n2
```

Compute **closest point between lines**.

Midpoint of closest segment approximates tube center.

Distance from point to center approximates `r`.

Procedure:

```
for many pairs (i,j):  
    compute closest point between normal lines  
    measure distance to each point  
    r\_estimate = average
```

Robust approach:

```
collect many r estimates  
take median
```

This gives a very good initial `r`.

---

# Step 4 — Compute tube center points

For each vertex:

```
k\_i = p\_i + r \* n\_i
```

These lie on the **major circle**.

---

# Step 5 — Fit circle to k points

Now we just need to fit a **3D circle**.

Algorithm:

1.  Compute PCA of `k_i`
    
2.  The **smallest eigenvector = circle normal = torus axis**
    

Then project onto the circle plane:

```
u,v = orthonormal basis perpendicular to axis
```

Project:

```
x\_i = dot(k\_i - mean\_k, u)  
y\_i = dot(k\_i - mean\_k, v)
```

Fit 2D circle:

```
(x - cx)^2 + (y - cy)^2 = R^2
```

Circle center in 3D:

```
C = mean\_k + cx\*u + cy\*v
```

Major radius:

```
R = fitted circle radius
```

---

# Step 6 — Output torus

Parameters:

```
center = C  
axis   = a  
major\_radius = R  
minor\_radius = r
```

---

# Step 7 — Optional nonlinear refinement

Now refine all parameters with LM.

Residual for point `p`:

```
d = p - C  
axial = dot(d,a)  
radial = norm(d - axial\*a)  
  
ri = sqrt((radial - R)^2 + axial^2) - r
```

Minimize:

```
Σ ri²
```

This usually converges in <20 iterations.

---

# Minimum number of points

Theoretical minimum:

```
7
```

But for this normal-based method you want:

```
~15–20 points
```

because radius estimation uses many pairs.

---

# Why this works well for fillets

Fillets typically:

-   cover small torus arcs
    
-   have accurate normals
    
-   have nearly constant minor radius
    

Direct torus fitting often fails because the arc is small.

This method works because:

```
normals → tube centers → major circle
```

It reconstructs the geometry in two simpler steps.

---

# Degenerate behavior

This method **naturally handles cylinders and spheres**:

### Cylinder

Normals intersect the **axis line** instead of a circle.

Result:

```
R → ∞
```

Detect by:

```
circle radius >> bounding box
```

Then treat as cylinder.

---

### Sphere

Normals intersect at a **single point**.

So all:

```
k\_i ≈ constant
```

Detect by:

```
variance(k\_i) very small
```

Then fit sphere.

---

# Complexity

```
normal intersection pairs: O(N)  
circle fit: O(N)  
LM refinement: O(kN)
```

Very fast.

---

# Practical improvements

Use:

```
random subset of pairs  
robust median radius  
reject pairs with nearly parallel normals
```

Also discard points where:

```
|dot(n, axis)| > 0.9
```

because those lie near the torus inner/outer extremes.

---

✅ **Summary**

Robust torus detection pipeline:

```
1\. PCA of normals → axis  
2\. intersect normal pairs → estimate r  
3\. k\_i = p\_i + r n\_i  
4\. circle fit of k\_i → center + R  
5\. optional LM refinement
```

This is the approach used in many **scan-to-CAD / reverse-engineering systems** because it handles **small toroidal patches** far better than direct torus fitting.

---

If you'd like, I can also show you a **very powerful extension used in industrial CAD reconstruction**: detecting toroidal fillets **directly from curvature tensors of the mesh**. It’s extremely effective for automatically finding fillets in large meshes.

You: 

yes

ChatGPT: 

Nice — this is where things start to look like the algorithms used in commercial **scan-to-CAD** systems. The key idea is that **tori have a very distinctive curvature signature**, so you can detect them *before* fitting anything heavy.

I'll give you a practical implementation-level description.

---

# 1\. Key geometric property of a torus

At each surface point there are two principal curvatures:

-   $k_1$
    
-   $k_2$
    

For a torus:

$$
k_1 = \frac{1}{r}
$$
 
$$
k_2 = \frac{\cos\theta}{R + r\cos\theta}
$$

Where:

-   $r$ = minor radius
    
-   $R$ = major radius
    
-   $\theta$ = angle around the tube
    

Important observations:

### First curvature

```
Plain text

k1 ≈ constant = 1/r
```

Across the entire torus.

### Second curvature

```
Plain text

k2 varies
```

It changes sign across the torus interior/exterior.

So toroidal patches look like:

| property | torus |
| --- | --- |
|  | k1 |
| k2 | varies smoothly |
| principal directions | rotate around axis |

That constant $k_1$ is the **fillet radius**.

---

# 2\. Compute curvature on a triangulated mesh

You need principal curvature estimates.

Standard method:

### Quadric surface fitting

For each vertex:

1.  gather neighbor vertices (1-ring or radius)
    
2.  build local tangent frame
    
3.  fit quadratic surface
    

Local surface:

```
z = ax² + bxy + cy²
```

From this you get:

```
k1, k2  
principal directions
```

Libraries implementing this:

-   libigl
    
-   CGAL
    
-   PCL
    

Or implement the common **Taubin / Rusinkiewicz curvature estimator**.

For implementation simplicity:

```
Plain text

Rusinkiewicz 2004 curvature estimator
```

is widely used for meshes.

---

# 3\. Detect toroidal candidate regions

For each vertex:

Compute

```
Plain text

k\_large = max(|k1|, |k2|)  
k\_small = min(|k1|, |k2|)
```

For toroidal fillets:

```
Plain text

k\_large ≈ constant  
k\_small varies
```

Define tolerance:

```
|k\_large - median(k\_large)| < ε
```

Typical ε:

```
ε ≈ 5–10%
```

Cluster connected vertices satisfying this.

Each cluster is a **candidate fillet patch**.

---

# 4\. Estimate minor radius

For each cluster:

```
r = 1 / median(k\_large)
```

This is extremely stable.

Much more stable than normal intersections.

---

# 5\. Estimate torus axis

Use **principal directions**.

The direction associated with **small curvature** points roughly along the torus major circle.

Procedure:

For each vertex:

```
t\_i = principal direction of smaller curvature
```

These directions lie approximately tangent to the major circle.

Compute the axis as the direction minimizing:

```
dot(axis, t\_i) ≈ 0
```

In practice:

```
axis = PCA smallest eigenvector of t\_i
```

Normalize.

---

# 6\. Compute tube center points

Same trick as earlier.

For each vertex:

```
k\_i = p\_i + r \* n\_i
```

These lie on the major circle.

---

# 7\. Fit major circle

Exactly the same circle fitting described earlier.

Steps:

```
1\. PCA of k\_i → circle plane  
2\. project onto plane  
3\. 2D circle fit
```

Results:

```
center C  
major radius R  
axis a
```

---

# 8\. Optional refinement

Run LM minimizing torus residual:

```
d = p - C  
axial = dot(d,a)  
radial = norm(d - axial\*a)  
  
ri = sqrt((radial - R)^2 + axial^2) - r
```

Usually only small adjustments needed.

---

# 9\. How to distinguish primitives

This curvature method also classifies surfaces automatically.

| surface | curvature signature |
| --- | --- |
| plane | k1≈0, k2≈0 |
| sphere | k1≈k2 constant |
| cylinder | k1 constant, k2≈0 |
| cone | k1 varies, k2≈0 |
| torus | k1 constant, k2 varies |

So you can build a primitive detector:

```
compute curvature  
classify surface type  
fit primitive
```

This is the architecture used in most **automatic feature recognition** systems.

---

# 10\. Minimum patch size

For reliable torus detection:

```
≥ 30–50 vertices
```

because curvature estimation is noisy.

But this is fine for fillets.

---

# 11\. Why this works extremely well for CAD meshes

CAD fillets have:

```
constant radius  
smooth normals  
clean triangulation
```

So:

```
k1 ≈ constant
```

is a very strong signal.

This makes toroidal fillets **easy to detect automatically**.

---

# 12\. Full torus detection pipeline

Putting everything together:

```
1\. compute vertex normals  
2\. estimate principal curvatures  
3\. classify surface regions  
4\. extract torus candidate clusters  
5\. r = 1 / median(|k\_large|)  
6\. axis = PCA of principal directions  
7\. k\_i = p\_i + r\*n\_i  
8\. fit circle to k\_i → C,R  
9\. optional nonlinear refinement
```

---

# 13\. Performance

Typical cost for mesh with N vertices:

```
curvature estimation: O(N)  
region clustering: O(N)  
circle fitting: O(M)  
LM refinement: O(kM)
```

Very fast.

---

# 14\. Practical tip for fillets

If you're detecting **fillets along edges**, you can get an even stronger signal:

The torus axis is approximately **parallel to the edge curve** of the two adjacent surfaces.

Some systems detect fillets by analyzing **dihedral angles along edges**, then fitting torus patches only there.

---

✅ If you'd like, I can also show you something extremely useful for your application:

**How commercial CAD kernels detect fillets in meshes in ~O(N)** using a clever trick involving **normal cones and medial axis estimation.**

That method is shockingly robust and avoids curvature estimation entirely.

You: 

yes

ChatGPT: 

This next trick is one of the **most robust fillet detectors** used in scan-to-CAD pipelines because it avoids computing curvature (which is noisy on meshes). Instead it uses **normal cones and medial geometry**. The idea is surprisingly simple once you see it.

I’ll describe a version that’s straightforward to implement.

---

# 1\. Key geometric property of a fillet torus

A fillet between two surfaces is a **constant-radius rolling surface**.

At every point:

```
Plain text

p + r n = tube center point
```

Where

-   `p` = surface point
    
-   `n` = surface normal
    
-   `r` = fillet radius
    

All those tube-center points lie on the **spine curve** (the torus major circle or arc).

So if you know the radius, you can map the surface to the spine directly.

The trick is **finding r without fitting a torus first**.

---

# 2\. Observation about normals on fillets

For two nearby vertices:

```
Plain text

p1, n1  
p2, n2
```

On a constant-radius surface, their normal lines intersect approximately at distance `r`.

Lines:

```
L1 = p1 + t n1  
L2 = p2 + s n2
```

Their closest approach gives the center of the osculating sphere.

Distance to points ≈ `r`.

This gives **local radius estimates**.

---

# 3\. Build a radius histogram

For many nearby vertex pairs:

Algorithm:

```
for vertex i  
    for neighbor j  
        compute closest points between normal lines  
        compute midpoint m  
        r\_estimate = |m - p\_i|  
        store r\_estimate
```

Reject pairs where normals are nearly parallel:

```
|dot(n\_i, n\_j)| > 0.95
```

After collecting many estimates:

```
r = median(r\_estimate)
```

If a strong peak exists → constant-radius surface detected.

This works extremely well for fillets.

---

# 4\. Compute tube-center points

Now compute:

```
k\_i = p\_i + r \* n\_i
```

These lie near the **fillet spine curve**.

For a toroidal fillet:

```
k\_i → arc of circle
```

For a cylindrical fillet:

```
k\_i → straight line
```

For a spherical blend:

```
k\_i → single point
```

This gives a **unified primitive detection**.

---

# 5\. Fit the spine geometry

Now simply fit geometry to `k_i`.

### If points form a line

Fit 3D line via PCA.

→ **cylindrical fillet**

### If points form a circle

Fit 3D circle.

→ **toroidal fillet**

### If points cluster to point

Compute centroid.

→ **spherical blend**

Detection metric:

```
λ1 ≥ λ2 ≥ λ3
```

Eigenvalues of covariance of `k_i`.

| shape | eigenvalue pattern |
| --- | --- |
| sphere | λ1≈λ2≈λ3≈0 |
| cylinder | λ1≫λ2≈λ3 |
| torus | λ1≈λ2≫λ3 |

---

# 6\. Recover torus parameters

If torus detected:

Circle fit on `k_i` gives:

```
center C  
axis a  
major radius R
```

Fillet radius already known:

```
minor radius = r
```

Done.

---

# 7\. Optional nonlinear refinement

Use LM with torus residual:

```
d = p - C  
axial = dot(d,a)  
radial = norm(d - axial\*a)  
  
ri = sqrt((radial - R)^2 + axial^2) - r
```

Minimize:

```
Σ ri²
```

Usually only tiny corrections needed.

---

# 8\. Why this works extremely well

Because fillets satisfy **two very strong constraints**:

1.  constant radius
    
2.  normals pointing toward spine
    

So the algorithm reduces the 3D surface problem to:

```
surface → medial curve
```

which is much simpler.

---

# 9\. Minimum points required

Technically:

```
~10 points
```

But stable detection:

```
30–50 mesh vertices
```

Typical fillet patches have hundreds.

---

# 10\. Complexity

For mesh with `N` vertices:

```
normal intersection estimates  O(N)  
spine fitting                  O(N)  
optional LM                    O(kN)
```

Very fast.

---

# 11\. Full fillet detection pipeline

```
1\. compute vertex normals  
2\. estimate radius from normal intersections  
3\. detect constant-radius clusters  
4\. compute k\_i = p\_i + r\*n\_i  
5\. classify k\_i geometry  
6\. fit spine (line/circle/point)  
7\. recover primitive parameters  
8\. optional nonlinear refinement
```

This pipeline detects:

-   edge fillets
    
-   rounded corners
    
-   cylindrical blends
    
-   spherical blends
    

using one mechanism.

---

# 12\. Practical robustness tricks

Use:

```
random subset of vertex pairs  
median radius  
RANSAC for circle fit
```

Also discard vertices where:

```
|dot(n\_i, average\_normal)| > 0.95
```

These are near extremal points where radius estimates become unstable.

---

✅ **Why this is powerful**

It transforms the problem:

```
fit torus surface
```

into

```
estimate medial curve of constant-radius surface
```

which is dramatically more stable.

This is very close to what tools like:

-   Geomagic
    
-   PolyWorks
    
-   Siemens NX Reverse Engineering
    

do internally.

---

If you're building a **mesh → CAD primitive detector**, I can also show you something extremely useful next:

**A single unified primitive-fitting framework (planes, spheres, cylinders, cones, tori) based on algebraic quadric fitting** that many geometry kernels use internally.

You: 

yes

ChatGPT: 

Nice—this gets into a powerful idea used in many geometry systems: **fit a general quadric first**, then **classify the surface type** (plane, sphere, cylinder, cone, etc.). It gives you a **single fitting pipeline** instead of separate solvers for each primitive.

Important caveat up front:

-   **Planes, spheres, cylinders, and cones are quadrics**
    
-   **Tori are *not* quadrics** (they’re quartic surfaces)
    

So the unified quadric fit works for **most primitives**, and torus detection typically happens **afterwards** (often using the constant-radius tricks we discussed).

---

# 1\. General quadric surface

A quadric in 3D is:

```
Ax² + By² + Cz²  
\+ Dxy + Exz + Fyz  
\+ Gx + Hy + Iz  
\+ J = 0
```

There are **10 coefficients**, but scale is arbitrary → **9 DOF**.

For a point `(x,y,z)` define the vector:

```
v = \[x², y², z², xy, xz, yz, x, y, z, 1\]
```

Then:

```
v · q = 0
```

where

```
q = \[A,B,C,D,E,F,G,H,I,J\]
```

---

# 2\. Least-squares quadric fitting

Given N points:

Build matrix:

```
M\_i = \[x\_i², y\_i², z\_i², x\_i y\_i, x\_i z\_i, y\_i z\_i, x\_i, y\_i, z\_i, 1\]
```

Solve:

```
M q ≈ 0
```

Minimize:

```
||M q||²
```

subject to:

```
||q|| = 1
```

Solution:

```
q = eigenvector of (MᵀM) with smallest eigenvalue
```

or simply:

```
SVD(M)
```

Take the right singular vector with smallest singular value.

Minimum points required:

```
9
```

but use ≥20 in practice.

---

# 3\. Convert quadric to matrix form

Represent the quadric as:

```
Q(x) = xᵀ A x + bᵀ x + c
```

Where:

```
A =  
\[ A  D/2 E/2  
  D/2 B  F/2  
  E/2 F/2 C \]  
  
b = \[G,H,I\]  
  
c = J
```

This representation is convenient for classification.

---

# 4\. Classify the quadric

Compute eigenvalues of matrix `A`.

The **rank and signs** determine the primitive.

### Plane

```
rank(A) = 0
```

Equation reduces to:

```
Gx + Hy + Iz + J = 0
```

---

### Sphere

Conditions:

```
A ≈ B ≈ C  
D ≈ E ≈ F ≈ 0
```

Meaning:

```
A ≈ λ I
```

Then center:

```
center = -½ A⁻¹ b
```

Radius:

```
r² = centerᵀ A center - c
```

---

### Cylinder

One eigenvalue of `A` ≈ 0.

Interpretation:

```
curvature in two directions  
flat in one
```

Axis = eigenvector corresponding to near-zero eigenvalue.

Radius computed from remaining terms.

---

### Cone

Matrix `A` full rank but determinant ≈ 0.

This indicates a quadratic surface passing through a vertex.

---

### Ellipsoid

All eigenvalues same sign.

---

### Hyperboloid / paraboloid

Mixed signs.

---

# 5\. Extract primitive parameters

Example: cylinder.

Steps:

1.  Eigen-decompose `A`.
    

```
A = R Λ Rᵀ
```

2.  Transform coordinates to principal frame.
    
3.  Reduce equation to canonical form:
    

```
x² + y² = r²
```

4.  Transform back to world coordinates.
    

---

# 6\. Residual error

For point `p` evaluate:

```
Q(p)
```

But this is **algebraic error**, not geometric distance.

Better metric:

```
|Q(p)| / ||∇Q(p)||
```

This approximates **point-to-surface distance**.

---

# 7\. Unified primitive detection pipeline

```
1\. collect mesh patch points  
2\. fit quadric via SVD  
3\. classify quadric type  
4\. extract parameters  
5\. compute geometric residual
```

If residual small → accept primitive.

Otherwise:

```
try torus detection
```

---

# 8\. Where torus fits in

Because torus is **quartic**, quadric fitting usually produces a **poor residual**.

So your pipeline becomes:

```
plane?  
sphere?  
cylinder?  
cone?  
quadric residual high → try torus detection
```

This avoids expensive torus fitting unless necessary.

---

# 9\. Advantages of quadric fitting

Single fast solve:

```
O(N)
```

No nonlinear optimization needed.

Also robust to noise.

Works very well for:

-   planes
    
-   spheres
    
-   cylinders
    
-   cones
    
-   paraboloids
    

---

# 10\. Minimum points

Theoretical:

```
9
```

Recommended:

```
20–50
```

---

# 11\. Numerical stability tips

Normalize coordinates first:

```
translate centroid to origin  
scale bounding box to ~1
```

Otherwise `x²` terms become very large.

After fitting, transform parameters back.

---

# 12\. Why most CAD reconstruction systems do this

Because it allows:

```
primitive detection = linear algebra
```

instead of nonlinear optimization.

Nonlinear solvers only used for:

-   torus
    
-   blends
    
-   freeform surfaces
    

---

✅ **Typical industrial primitive detection pipeline**

```
1\. cluster mesh patches  
2\. fit quadric  
3\. classify primitive  
4\. if not quadric → try torus  
5\. if not torus → fit NURBS
```

This is extremely fast and scales to large meshes.

---

If you're interested, I can also show you something **very relevant to your OCCT work**:

There is a **beautiful trick to fit cylinders, cones, and spheres using only normals (no positions)** that often produces dramatically more stable fits on triangulated CAD meshes.

You: 

yes

ChatGPT: 

This is a surprisingly elegant trick: **many primitives can be identified almost entirely from normals**, because their normals obey very simple geometric constraints. For meshes derived from CAD, normals are usually very accurate, so these fits can be **more stable than point-based fits**.

I'll show the practical algorithms for:

-   plane
    
-   sphere
    
-   cylinder
    
-   cone
    

All using normals.

---

# 1\. Plane (trivial)

For a plane:

```
normal is constant everywhere
```

Algorithm:

```
n = normalize(mean(normals))
```

Plane offset:

```
d = -dot(n, centroid(points))
```

Equation:

```
n·x + d = 0
```

Minimum points:

```
3
```

---

# 2\. Sphere (normal intersection)

Sphere property:

```
surface normals intersect at sphere center
```

For each vertex:

```
p\_i + t n\_i
```

should pass through center `C`.

Solve:

```
p\_i + t\_i n\_i = C
```

Rearrange:

```
C - p\_i = t\_i n\_i
```

Take cross product with `n_i`:

```
(C - p\_i) × n\_i = 0
```

Which expands to linear equations in `C`.

For each point:

```
\[  0  -n\_z  n\_y \]   \[Cx - px\]  
\[ n\_z   0  -n\_x \] \* \[Cy - py\] = 0  
\[-n\_y  n\_x   0  \]   \[Cz - pz\]
```

Stack for many points and solve least squares.

Then:

```
radius = mean(|p\_i - C|)
```

Minimum points:

```
4
```

---

# 3\. Cylinder (normal cone trick)

Key property:

```
cylinder normals are perpendicular to axis
```

So for axis `a`:

```
dot(n\_i, a) = 0
```

Stack equations:

```
N a = 0
```

Where:

```
N =  
\[n1x n1y n1z  
 n2x n2y n2z  
 ...  
\]
```

Solve with SVD.

Axis = singular vector with smallest singular value.

---

## Cylinder center

After axis known:

For each point:

```
d\_i = p\_i - (p\_i·a)a
```

Project into plane perpendicular to axis.

Now it's a **2D circle fit**.

Steps:

```
1\. choose basis u,v ⟂ a  
2\. project points  
3\. fit circle
```

This gives:

```
cylinder axis point  
radius
```

Minimum points:

```
5–6
```

---

# 4\. Cone (normals meet on sphere)

Cone property:

```
normals intersect along a cone apex sphere
```

Another useful property:

```
all normals make constant angle with axis
```

So:

```
dot(n\_i, a) = cos(θ)
```

Algorithm:

### Step 1: find axis

Compute PCA of normals.

Axis ≈ mean normal direction.

### Step 2: estimate cone angle

```
θ = acos(mean(dot(n\_i, a)))
```

### Step 3: find apex

Each point satisfies:

```
(p\_i - A) · n\_i = 0
```

where `A` = apex.

This forms linear system:

```
n\_i · A = n\_i · p\_i
```

Solve least squares.

Minimum points:

```
6–8
```

---

# 5\. Why normal-based fitting is powerful

Normals remove one dimension of noise.

Point-based fitting uses:

```
position errors
```

Normal-based fitting uses:

```
orientation constraints
```

Which are often **more stable on triangulated CAD meshes**.

Example:

| primitive | best fitting method |
| --- | --- |
| plane | normals |
| sphere | normals |
| cylinder | normals + points |
| cone | normals |
| torus | normals + medial method |

---

# 6\. Combined primitive detection pipeline

This is a very effective architecture:

```
1\. compute vertex normals  
2\. cluster mesh patches  
3\. detect primitive type from normals  
4\. compute parameters using normals  
5\. refine with point-based least squares
```

This hybrid approach gives:

```
fast detection  
accurate parameters
```

---

# 7\. Degenerate detection

Normal distributions also classify primitives:

### Plane

```
normals cluster tightly
```

### Sphere

```
normals distributed over sphere
```

### Cylinder

```
normals lie on great circle
```

### Cone

```
normals lie on small circle
```

These patterns are extremely useful for classification.

---

# 8\. Minimum points summary

| primitive | minimum |
| --- | --- |
| plane | 3 |
| sphere | 4 |
| cylinder | 5–6 |
| cone | 6–8 |
| torus | ~20 |

---

# 9\. Practical tip for CAD meshes

Before fitting primitives, do **region growing by normal similarity**:

```
|dot(n\_i, n\_j)| > 0.95
```

This gives clean patches that correspond to analytic surfaces.

Then run primitive detection per patch.

This dramatically improves robustness.

---

✅ If you're building a **mesh → analytic surface detector**, the next thing that would probably help you most is a **complete architecture used in CAD kernels for reconstructing analytic surfaces from triangulations** (including patch segmentation, primitive ranking, and fallback to NURBS). It's a really nice pipeline and fits perfectly with the algorithms we've been discussing.

You: 

yes

ChatGPT: 

What follows is a **practical architecture for mesh → analytic surface reconstruction** that’s very close to what reverse-engineering systems and CAD kernels use internally. It’s designed to be:

-   **robust to noisy triangulations**
    
-   **fast (mostly linear algebra)**
    
-   **modular (add primitives easily)**
    

I’ll describe the pipeline in a way you could realistically implement.

---

# Overview

Goal: convert a triangulated mesh into analytic patches:

```
plane  
sphere  
cylinder  
cone  
torus  
NURBS fallback
```

Pipeline:

```
1\. compute normals + curvature  
2\. segment mesh into candidate patches  
3\. classify primitive type  
4\. estimate primitive parameters  
5\. refine parameters  
6\. validate residual error  
7\. grow patch  
8\. if fitting fails → NURBS
```

---

# 1\. Preprocessing

Compute per-vertex:

```
normal  
principal curvature (optional)
```

Normals:

```
area-weighted face average
```

Curvature:

Rusinkiewicz estimator works well.

This information drives segmentation.

---

# 2\. Mesh segmentation

Goal: break mesh into **smooth regions**.

Region growing algorithm:

Seed triangle.

Add neighbor if:

```
|dot(n\_i, n\_seed)| > cos(angle\_threshold)
```

Typical threshold:

```
10–15°
```

Also require curvature similarity if available.

Result:

```
patches of roughly constant curvature
```

These correspond roughly to CAD faces.

---

# 3\. Primitive classification

For each patch, analyze **normal distribution**.

Compute covariance of normals:

```
C = Σ n\_i n\_iᵀ
```

Eigenvalues:

```
λ1 ≥ λ2 ≥ λ3
```

Patterns:

| primitive | normal pattern |
| --- | --- |
| plane | λ1≫λ2≈λ3 |
| cylinder | λ1≈λ2≫λ3 |
| sphere | λ1≈λ2≈λ3 |
| cone | λ1≈λ2≫λ3 but normals tilted |
| torus | more complex |

Also use curvature:

```
k1,k2
```

Rules:

```
if k1≈k2 ≠0 → sphere  
if k1≈const and k2≈0 → cylinder  
if k1≈const and k2 varies → torus
```

This step narrows candidate primitives.

---

# 4\. Parameter estimation

Use specialized solvers.

### Plane

```
PCA of points
```

### Sphere

```
normal intersection
```

### Cylinder

```
axis from normals  
circle fit
```

### Cone

```
constant-angle normals  
solve apex
```

### Torus

```
estimate minor radius  
compute k\_i = p\_i + r n\_i  
fit circle to k\_i
```

These produce initial parameters.

---

# 5\. Nonlinear refinement

Now refine parameters using geometric residuals.

Example: cylinder

```
ri = distance(point, axis) − r
```

Example: sphere

```
ri = |p − C| − r
```

Example: torus

```
ri = sqrt((radial − R)^2 + axial^2) − r
```

Minimize:

```
Σ ri²
```

Use:

```
Levenberg–Marquardt
```

Usually <20 iterations.

---

# 6\. Residual validation

Compute:

```
RMS error  
max error
```

Accept primitive if:

```
RMS < tolerance  
max < 3× tolerance
```

Tolerance typically:

```
0.1–0.5% of bounding box
```

---

# 7\. Patch growing

Once primitive found, expand patch.

For neighboring triangles:

Add if point distance to primitive < tolerance.

This grows the analytic face beyond the original seed cluster.

---

# 8\. Boundary extraction

Find boundary edges between patches.

These correspond to CAD **edge curves**.

Fit curves:

```
line  
circle  
ellipse  
BSpline
```

---

# 9\. Topology reconstruction

Build CAD topology:

```
faces → surfaces  
edges → curves  
vertices → curve intersections
```

Ensure adjacency consistency.

For example:

```
two cylinders intersect → circle edge  
plane-cylinder → line edge
```

---

# 10\. NURBS fallback

If primitive fitting fails:

Fit NURBS surface.

Typical approach:

```
1\. parameterize patch  
2\. fit bicubic B-spline  
3\. adaptive knot insertion
```

Only needed for freeform areas.

---

# 11\. Final CAD model

Output:

```
analytic surfaces  
trim curves  
topological adjacency
```

For OCCT this becomes:

```
Geom\_Plane  
Geom\_CylindricalSurface  
Geom\_SphericalSurface  
Geom\_ToroidalSurface  
Geom\_ConicalSurface  
Geom\_BSplineSurface
```

Then assemble with:

```
TopoDS\_Face  
TopoDS\_Edge  
TopoDS\_Vertex
```

---

# 12\. Performance

For mesh with N vertices:

```
normal + curvature   O(N)  
segmentation         O(N)  
primitive fitting    O(N\_patch)  
LM refinement        small
```

Overall roughly:

```
O(N)
```

Which is why these pipelines scale to millions of triangles.

---

# 13\. Practical heuristics used in industry

A few tricks that improve results dramatically:

### ignore tiny patches

```
< 20 triangles
```

### prefer simpler primitives

Ranking:

```
plane  
sphere  
cylinder  
cone  
torus  
NURBS
```

### merge adjacent patches

If parameters nearly identical.

---

# 14\. Result quality

Good meshes typically reconstruct with:

```
\>90% analytic surfaces
```

Freeform areas become NURBS.

This matches how most CAD models are structured.

---

✅ Since you're working with **OCCT geometry and triangulated surfaces**, one more thing would likely help a lot:

There is a **very elegant algorithm used in CAD kernels to detect cylinders and cones directly from triangle meshes using Gaussian maps** (mapping normals to the unit sphere). It makes primitive classification almost trivial and extremely fast.

You: 

yes

ChatGPT: 

The **Gaussian map trick** is one of the cleanest ways to classify analytic surfaces from a triangle mesh. It’s widely used in geometry processing because it converts a messy 3D problem into a **simple shape recognition problem on the unit sphere**.

The idea:  
Map every surface normal to the **unit sphere**. The pattern formed by those points reveals the surface type immediately.

---

# 1\. Gaussian map

For each mesh vertex or triangle:

```
Plain text

n = surface normal (unit)
```

Plot `n` as a point on the **unit sphere**.

This is the **Gaussian map**.

So a mesh patch becomes a **point cloud on S²**.

Different analytic surfaces produce distinctive patterns.

---

# 2\. Patterns for common primitives

### Plane

All normals identical.

On sphere:

```
Plain text

single point
```

Detection:

```
variance(normals) ≈ 0
```

---

### Sphere

Normals cover a full spherical region.

On sphere:

```
Plain text

patch on sphere
```

Detection:

```
normals distributed roughly uniformly
```

Also:

```
n ≈ normalize(p − center)
```

---

### Cylinder

Cylinder normals rotate around the axis.

Key property:

```
dot(n, axis) = 0
```

On the Gaussian sphere:

```
Plain text

great circle
```

Because normals lie in plane through origin.

Detection:

Fit plane through normal points:

```
axis = plane normal
```

---

### Cone

Cone normals satisfy:

```
dot(n, axis) = cos(θ)
```

On Gaussian sphere:

```
Plain text

small circle
```

Detection:

Fit circle on sphere.

---

### Torus

Torus normals form a **band** between two circles.

On sphere:

```
Plain text

annular region
```

More complex, but still structured.

---

# 3\. Practical detection algorithm

Given normals `n_i`.

Compute covariance:

```
C = Σ n\_i n\_iᵀ
```

Eigenvalues:

```
λ1 ≥ λ2 ≥ λ3
```

Interpretation:

| surface | eigenvalues |
| --- | --- |
| plane | λ1≫λ2≈λ3≈0 |
| cylinder | λ1≈λ2≫λ3 |
| sphere | λ1≈λ2≈λ3 |
| cone | λ1≈λ2≫λ3 but shifted |

This alone already classifies many patches.

---

# 4\. Explicit great-circle fitting (cylinder detection)

Normals lie on plane:

```
a · n = 0
```

Where `a` = cylinder axis.

Solve:

```
N a = 0
```

with

```
N = matrix of normals
```

Use SVD.

Smallest singular vector = axis.

Very stable.

---

# 5\. Detecting cones via small circles

For cone normals:

```
a · n = cos θ
```

Algorithm:

1.  find axis via PCA
    
2.  compute
    

```
cosθ\_i = dot(a, n\_i)
```

3.  check variance.
    

If small → cone.

Angle:

```
θ = acos(mean(cosθ\_i))
```

---

# 6\. Detecting spheres from Gaussian map

Sphere property:

```
n ≈ normalize(p − C)
```

So normals intersect at center.

Compute center via normal intersection (as described earlier).

Verify:

```
|p − C| ≈ constant
```

---

# 7\. Detecting torus from Gaussian map

Torus normals satisfy:

```
dot(n, axis) ∈ \[−sin(φ), sin(φ)\]
```

So on sphere:

```
band around equator
```

Detection:

1.  estimate axis via PCA
    
2.  compute
    

```
d\_i = dot(n\_i, axis)
```

3.  check if values span symmetric interval.
    

Then run torus fit.

---

# 8\. Why Gaussian maps are powerful

They remove translation and scale.

Normals encode only **surface orientation**, which depends only on surface type.

So classification becomes:

```
pattern recognition on sphere
```

instead of 3D fitting.

This is extremely robust.

---

# 9\. Full primitive classification using Gaussian map

Algorithm:

```
1\. compute normals  
2\. map normals to unit sphere  
3\. compute covariance of normals  
4\. inspect eigenvalues
```

Decision tree:

```
if variance ≈0 → plane  
else if λ1≈λ2≈λ3 → sphere  
else if λ1≈λ2≫λ3 → cylinder or cone  
else → torus or freeform
```

Then run specialized parameter solvers.

---

# 10\. Minimum points

Reliable classification usually needs:

```
20–30 normals
```

Which corresponds to roughly 20–30 mesh triangles.

---

# 11\. Why CAD meshes work well

Meshes exported from CAD usually have:

```
exact normals  
clean curvature  
regular triangulation
```

So Gaussian-map classification becomes extremely reliable.

---

# 12\. Integration into your primitive pipeline

Combine with earlier steps:

```
segmentation  
↓  
Gaussian-map classification  
↓  
primitive-specific solver  
↓  
LM refinement
```

This reduces expensive fitting dramatically.

---

# 13\. Extra trick: normal-space clustering

Sometimes segmentation can be done **entirely in normal space**.

Cluster normals using:

```
k-means on sphere
```

Each cluster corresponds to a surface patch.

This is surprisingly effective for CAD meshes.

---

✅ If you're building a **high-quality mesh → CAD reconstruction tool**, the final trick that would probably help the most is something used in industrial kernels:

**RANSAC-based primitive detection combined with Gaussian-map filtering**, which allows detecting primitives even in messy meshes or partial scans. It’s a very elegant hybrid approach.

You: 

Forget about normal-based methods for the moment - describe how to implement conical fitting to a set of points, similar to the earlier algorithms for toroids and cylinders.

ChatGPT: 

Here is a **practical cone-fitting algorithm from points only**, analogous to the cylinder and torus fits we discussed earlier. The approach is:

1.  estimate axis
    
2.  estimate apex
    
3.  estimate cone angle
    
4.  refine with nonlinear least squares
    

This is stable and relatively simple to implement.

---

# 1\. Cone parameterization

A cone is defined by:

```
apex A : Vec3  
axis a : unit Vec3  
half-angle θ
```

For any point `p` on the cone:

```
angle between (p − A) and a = θ
```

Equivalently:

```
dot( normalize(p − A), a ) = cos θ
```

For fitting we avoid normalization.

Define:

```
v = p − A
```

Residual:

```
ri = dot(v, a)^2 − |v|^2 cos²θ
```

This equals zero for perfect points.

We minimize:

```
Σ ri²
```

---

# 2\. Minimum number of points

Unknown parameters:

```
A (3)  
axis direction (2 dof)  
θ (1)
```

Total:

```
6 parameters
```

Minimum points:

```
6
```

But stable fits require:

```
20–30 points
```

---

# 3\. Step 1 — estimate axis

Compute PCA of points.

Let eigenvalues:

```
λ1 ≥ λ2 ≥ λ3
```

For a cone:

```
axis ≈ eigenvector with smallest eigenvalue
```

Reason: points spread mainly around the axis.

Normalize:

```
a = normalize(eigenvector)
```

---

# 4\. Step 2 — estimate apex

The apex lies on the axis.

Let:

```
A = c + t a
```

Where:

```
c = centroid of points
```

So only unknown is `t`.

For each point:

```
v\_i = p\_i − (c + t a)
```

Cone condition:

```
dot(v\_i, a)^2 = |v\_i|^2 cos²θ
```

But we don't know θ yet.

Instead compute **distance from axis**.

Define:

```
h\_i = dot(p\_i − A, a)  
r\_i = |(p\_i − A) − h\_i a|
```

For a cone:

```
r\_i = h\_i tan θ
```

This is linear if A known.

But since `A = c + t a`, solve for `t` that best satisfies:

```
r\_i / h\_i ≈ constant
```

Practical method:

Search `t` along axis that minimizes variance of:

```
k\_i = r\_i / h\_i
```

Then:

```
tan θ = median(k\_i)
```

So:

```
θ = atan(median(k\_i))
```

---

# 5\. Step 3 — compute apex

Once `t` chosen:

```
A = c + t a
```

Now compute:

```
h\_i = dot(p\_i − A, a)  
r\_i = |(p\_i − A) − h\_i a|
```

Angle:

```
θ = atan( mean(r\_i / h\_i) )
```

---

# 6\. Step 4 — nonlinear refinement

Now refine parameters:

```
A  
a  
θ
```

Use LM.

Residual:

```
v = p − A  
ri = dot(v,a)^2 − |v|^2 cos²θ
```

Alternative (more geometric):

```
ri = angle(v,a) − θ
```

But angle form is slower.

---

# 7\. Implementation outline

Pseudo-pipeline:

```
1\. compute centroid  
2\. PCA → axis estimate  
3\. search t along axis to estimate apex  
4\. compute θ  
5\. LM refinement
```

---

# 8\. Detecting degenerate cases

Cone fitting often collapses into other primitives.

### Cylinder

When:

```
θ ≈ 90°
```

### Plane

When:

```
θ ≈ 0
```

### Sphere

When apex lies near centroid and angle varies.

So after fitting check:

```
θ small → plane  
θ near π/2 → cylinder
```

---

# 9\. Residual error

For point `p`:

```
v = p − A  
h = dot(v,a)  
r = |v − h a|
```

Expected:

```
r = h tan θ
```

Residual:

```
ri = r − h tan θ
```

Minimize:

```
Σ ri²
```

This is often the best residual.

---

# 10\. Complexity

```
PCA                O(N)  
axis search        O(N \* samples)  
LM refinement      O(kN)
```

Typically extremely fast.

---

# 11\. Numerical stability tips

Reject points where:

```
|h| small
```

because near the apex the ratio becomes unstable.

Also normalize coordinates to roughly unit scale before solving.

---

# 12\. Output

Final parameters:

```
struct ConeFit {  
    apex: Vec3  
    axis: Vec3  
    angle: f64  
    rms\_error: f64  
}
```

---

# 13\. Relation to other primitives

This structure fits nicely with earlier algorithms:

| primitive | parameters |
| --- | --- |
| cylinder | axis + radius |
| sphere | center + radius |
| cone | apex + axis + angle |
| torus | center + axis + R + r |

So your primitive detection framework can reuse the same **LM infrastructure**.

---

If you'd like, I can also show a **much better cone initialization method used in geometric modeling papers** that avoids the axis search entirely and solves the cone parameters **almost directly via linear algebra**. It’s surprisingly elegant and often produces a near-perfect initial fit.

You: 

yes

ChatGPT: 

Here’s a **much cleaner initialization method** for cone fitting that avoids searching along the axis. It comes from classic geometric modeling work (e.g., Eberly-style primitive fits). The key idea is to use the **implicit quadratic form of a cone** and solve most of the parameters with linear algebra before any nonlinear refinement.

This produces a **very good initial estimate** that LM can refine quickly.

---

# 1\. Implicit equation of a cone

A right circular cone with:

-   apex $A$
    
-   axis $a$ (unit)
    
-   half-angle $θ$
    

satisfies:

$$
((p-A)\cdot a)^2 = |p-A|^2 \cos^2 θ
$$

Let:

```
v = p - A
```

Then:

```
(v·a)² − |v|² cos²θ = 0
```

Rewrite:

```
vᵀ (a aᵀ − cos²θ I) v = 0
```

Define matrix:

```
M = a aᵀ − cos²θ I
```

So the cone equation becomes:

```
(p − A)ᵀ M (p − A) = 0
```

This is a **quadric centered at the apex**.

---

# 2\. Expand the equation

Expand:

```
pᵀ M p − 2Aᵀ M p + Aᵀ M A = 0
```

Define unknowns:

```
M (symmetric 3×3)  
b = −2 M A  
c = Aᵀ M A
```

So:

```
pᵀ M p + bᵀ p + c = 0
```

This is a **general quadric equation**, but with the special constraint that it represents a cone.

---

# 3\. Linear least-squares solve

For point $p=(x,y,z)$:

Construct row:

```
\[x², y², z², xy, xz, yz, x, y, z, 1\]
```

Solve:

```
M q ≈ 0
```

as in quadric fitting via SVD.

This gives the coefficients:

```
\[A,B,C,D,E,F,G,H,I,J\]
```

representing:

```
Ax² + By² + Cz²  
\+ Dxy + Exz + Fyz  
\+ Gx + Hy + Iz  
\+ J = 0
```

---

# 4\. Extract cone parameters

Convert to matrix form:

```
Q(x) = xᵀ A x + bᵀ x + c
```

Where

```
A =  
\[ A  D/2 E/2  
  D/2 B  F/2  
  E/2 F/2 C \]  
  
b = \[G,H,I\]  
c = J
```

---

## Apex

For a cone the apex satisfies:

```
∇Q(A) = 0
```

Gradient:

```
2A A + b = 0
```

Solve:

```
A = −½ A⁻¹ b
```

This gives the **apex directly**.

---

## Axis direction

Compute eigen-decomposition of matrix `A`.

For a cone:

```
eigenvalues ≈ {λ, λ, μ}
```

The eigenvector corresponding to **distinct eigenvalue** gives the axis.

Normalize:

```
a = eigenvector
```

---

## Cone angle

Using relation:

```
M = a aᵀ − cos²θ I
```

The eigenvalues of `M` are:

```
{1 − cos²θ, −cos²θ, −cos²θ}
```

So:

```
cos²θ = −λ\_perp
```

Where `λ_perp` is one of the repeated eigenvalues.

Compute:

```
θ = acos( sqrt(cos²θ) )
```

---

# 5\. Result

You now have:

```
apex A  
axis a  
angle θ
```

This is usually already a good fit.

---

# 6\. Optional LM refinement

Refine parameters minimizing:

```
v = p − A  
h = dot(v,a)  
r = |v − h a|
```

Residual:

```
ri = r − h tan θ
```

Minimize:

```
Σ ri²
```

LM typically converges in ~10 iterations.

---

# 7\. Minimum number of points

The quadric solve has 9 DOF.

Minimum:

```
9 points
```

But use:

```
≥ 20
```

for stability.

---

# 8\. Degenerate cases

Check eigenvalues of matrix `A`.

| pattern | interpretation |
| --- | --- |
| λ≈λ≈0 | plane |
| λ≈λ≈λ | sphere |
| λ≈λ≫μ | cylinder |
| λ≈λ≠μ | cone |

This lets you **classify primitives automatically**.

---

# 9\. Implementation pipeline

```
1\. normalize coordinates  
2\. build quadric matrix  
3\. SVD → quadric coefficients  
4\. compute matrix A, vector b  
5\. solve apex = −½ A⁻¹ b  
6\. eigen-decompose A → axis  
7\. compute cone angle  
8\. LM refine
```

---

# 10\. Why this method is good

Advantages:

-   no axis search
    
-   purely linear initialization
    
-   very fast
    
-   works even for partial cones
    

This is why variants of this method appear in many **geometric modeling libraries**.

---

💡 Since you're already thinking about cylinders, cones, and tori, there’s one more algorithm that’s incredibly useful when building a primitive-fitting system:

A **single unified nonlinear solver that fits any surface of revolution** (cylinders, cones, spheres, tori) by solving for the **axis and a radial profile curve** simultaneously.

That approach can dramatically simplify the architecture of a reconstruction system.

You: 

yes

ChatGPT: 

This is a really nice trick used in some reverse-engineering / mesh-to-CAD pipelines: **fit all surfaces of revolution with the same solver**. Cylinders, cones, spheres, and tori all fall out as special cases.

The idea is:

> Instead of fitting the surface directly in 3D, first fit the **axis of revolution**, then reduce the problem to **2D curve fitting in the radial–axial plane**.

Once you see it, the implementation becomes surprisingly simple.

---

# 1\. Surfaces of revolution formulation

A surface of revolution is generated by rotating a **profile curve** around an axis.

Let

```
axis point  q  
axis dir    a (unit)
```

For any 3D point `p`, compute coordinates relative to the axis:

```
v = p - q  
h = dot(v, a)              // axial coordinate  
r = |v - h a|              // radial distance
```

So every 3D point becomes a **2D point**:

```
(h, r)
```

Now fitting a surface of revolution becomes:

```
fit curve r = f(h)
```

in the **(h,r) plane**.

---

# 2\. Profiles for common primitives

Each primitive corresponds to a simple curve in the `(h,r)` plane.

### Cylinder

```
r = R
```

horizontal line.

---

### Cone

```
r = (h - h0) tan θ
```

straight line.

---

### Sphere

If center is `(0,0)` in `(h,r)` coordinates:

```
r² + h² = R²
```

circle centered on axis.

---

### Torus

Profile is a circle offset from axis:

```
(r - R)² + h² = r\_minor²
```

circle whose center is at `(h=0, r=R)`.

---

# 3\. The unified algorithm

Pipeline:

```
1 estimate axis  
2 compute (h,r) coordinates  
3 fit candidate curves  
4 choose best primitive  
5 optional nonlinear refinement
```

This allows **one axis estimation + multiple cheap 2D fits**.

---

# 4\. Step 1 — axis estimation

Use PCA of points:

```
C = covariance(points)
```

Eigenvalues:

```
λ1 ≥ λ2 ≥ λ3
```

For surfaces of revolution:

```
axis ≈ eigenvector with smallest eigenvalue
```

Compute centroid:

```
q = mean(points)
```

Axis:

```
a = normalize(eigenvector\_min)
```

---

# 5\. Step 2 — compute reduced coordinates

For each point:

```
v = p - q  
h = dot(v,a)  
r = norm(v - h a)
```

Store:

```
(h\_i, r\_i)
```

Now you have a **2D dataset**.

---

# 6\. Fit primitives in 2D

## Cylinder

Fit:

```
r = R
```

Estimate:

```
R = mean(r\_i)
```

Residual:

```
ri = r\_i - R
```

---

## Cone

Fit line:

```
r = m h + b
```

Use standard linear regression.

Angle:

```
θ = atan(m)
```

Apex:

```
h\_apex = -b/m
```

---

## Sphere

Fit circle centered on axis:

```
h² + r² + Ah + B = 0
```

Solve linear least squares.

Center:

```
h0 = -A/2
```

Radius:

```
R = sqrt(h0² - B)
```

---

## Torus

Fit circle in `(h,r)` plane:

```
(h - hc)² + (r - Rc)² = r\_minor²
```

Use standard **circle least squares**.

Then:

```
major radius = Rc  
minor radius = r\_minor  
center shift = hc
```

---

# 7\. Choose best primitive

Compute RMS residual for each fit:

```
rms = sqrt(mean(ri²))
```

Choose the smallest.

Example:

```
cylinder rms = 0.002  
cone rms     = 0.0015  
sphere rms   = 0.02  
torus rms    = 0.0008
```

→ torus.

---

# 8\. Recover 3D parameters

Given `(h,r)` parameters.

### Torus

Center:

```
C = q + hc a
```

Axis:

```
a
```

Major radius:

```
R
```

Minor radius:

```
r
```

---

### Cone

Apex:

```
A = q + h\_apex a
```

Angle:

```
θ
```

---

### Sphere

Center:

```
C = q + h0 a
```

---

# 9\. Nonlinear refinement

Once primitive chosen, run a full 3D LM refinement.

Examples:

Cylinder:

```
ri = |(p - q) × a| - R
```

Cone:

```
ri = r\_i - h\_i tanθ
```

Torus:

```
ri = (sqrt((r\_i - R)² + h\_i²) - r\_minor)
```

---

# 10\. Minimum number of points

| primitive | minimum |
| --- | --- |
| cylinder | 5 |
| cone | 6 |
| sphere | 4 |
| torus | 7 |

But practical:

```
20+
```

---

# 11\. Why this method is extremely useful

Advantages:

**1\. Same pipeline for many primitives**

```
axis estimate  
→ 2D reduction  
→ cheap curve fits
```

---

**2\. Extremely fast**

All primitive detection becomes:

```
O(N)
```

after axis estimation.

---

**3\. Numerically stable**

Fitting circles and lines in 2D is far easier than solving full 3D nonlinear problems.

---

# 12\. Bonus: detecting fillets (tori) on CAD meshes

Fillet patches typically produce `(h,r)` points that lie on:

```
circle centered at r ≈ edge radius
```

Which makes torus detection in this framework **very reliable**.

---

💡 If you're building a **mesh → analytic surface reconstruction system**, the next step that dramatically improves robustness is something many people miss:

**Estimating the axis of revolution directly from the points without PCA**, using a method based on **intersection of normal planes**.

That technique is much more stable for **small torus patches and fillets**, where PCA often fails.

You: 

Given a large mesh, how to pick a small neighborhood of points to estimate an axis from?

ChatGPT: 

When you have a **large mesh**, the key problem is not the fitting algorithm — it's **choosing a patch of points that belong to the same surface**. Good axis estimation depends heavily on that patch being locally consistent.

The typical pipeline used in reverse-engineering / mesh-to-CAD tools is:

```
seed triangle  
→ grow local patch  
→ estimate primitive axis  
→ verify  
→ expand region
```

Below are practical methods that work well.

---

# 1\. Start from a seed triangle

Pick a triangle `t₀`.

Compute:

```
p₀ = triangle centroid  
n₀ = triangle normal
```

The seed triangle defines the **local orientation**.

---

# 2\. Build a local neighborhood

Use **mesh adjacency** to grow a small patch.

Breadth-first expansion:

```
patch = {t₀}  
  
while patch.size < N\_target:  
    add adjacent triangles
```

Typical sizes:

```
N\_target = 20–50 triangles
```

This corresponds to ~20–50 normals/points.

---

# 3\. Filter by normal similarity

Avoid crossing sharp edges.

For each candidate triangle:

```
if dot(n\_i, n₀) > cos(angle\_threshold)  
    accept
```

Typical:

```
angle\_threshold = 20–30°
```

This keeps the patch on one smooth surface.

---

# 4\. Filter by curvature consistency

For surfaces of revolution, curvature directions behave predictably.

Compute normal variation:

```
Δn = angle(n\_i, n₀)
```

Reject triangles where:

```
Δn grows too quickly with distance
```

This prevents jumping across corners.

---

# 5\. Gather points

Collect:

```
triangle centroids  
or  
triangle vertices
```

Centroids often work best.

So now you have:

```
points = {p\_i}  
normals = {n\_i}
```

---

# 6\. Estimate axis

Use the **normal-plane intersection method**.

For surfaces of revolution:

```
axis lies in every normal plane
```

Normal plane equation:

```
n\_i · (x − p\_i) = 0
```

Solve least squares for axis direction.

A practical approach:

Build matrix:

```
A\_i = I − n\_i n\_iᵀ
```

Solve:

```
min Σ |A\_i (x − p\_i)|²
```

The direction of minimal variance gives the axis.

This is **more stable than PCA**.

---

# 7\. Validate the axis

Once axis `(q,a)` estimated:

Compute `(h,r)` for each point.

```
v = p − q  
h = dot(v,a)  
r = |v − h a|
```

Check consistency.

Examples:

Cylinder:

```
var(r) small
```

Cone:

```
r / h approximately constant
```

Torus:

```
points lie on circle in (h,r)
```

If residual large → patch is not a surface of revolution.

---

# 8\. Expand region

If the fit is good, grow the region.

Add neighboring triangles if they satisfy:

```
distance\_to\_surface < tolerance
```

Example cylinder test:

```
| |(p − q) × a| − R | < tol
```

This becomes a **region growing primitive detector**.

---

# 9\. Stop conditions

Stop growing when:

```
residual exceeds tolerance  
or  
sharp edge encountered
```

You now have a full primitive patch.

---

# 10\. Repeat

Mark triangles as processed.

Pick another seed and repeat.

This is essentially how many **mesh segmentation → primitive fitting** systems work.

---

# 11\. Practical parameters

Good starting values:

```
seed patch size      30 triangles  
normal threshold     25°  
fit tolerance        0.1–1% of mesh size
```

---

# 12\. Performance

For mesh with `T` triangles:

```
axis estimation      O(N)  
patch growth         O(N)
```

Total detection often ends up roughly:

```
O(T)
```

---

# 13\. Important trick used in CAD reconstruction

Instead of picking **random seeds**, pick seeds where curvature is strong.

Compute per-vertex curvature estimate:

```
mean normal deviation
```

Seed at **high curvature areas** first.

This detects:

```
fillets  
cones  
spheres  
tori
```

very quickly.

---

💡 One trick that becomes extremely valuable when detecting **toroidal fillets around edges** (which sounds exactly like your use case) is to seed from **mesh edges where the normals rotate strongly around a common axis**.

That allows you to detect fillets with just **6–10 triangles**, instead of needing a larger patch.