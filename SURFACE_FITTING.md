# Surface Fitting Algorithms

This document describes the surface fitting math and validation algorithms used in brepper's Stage 2 pipeline. Each section covers the fitting math and validation for one surface type.

For the overall pipeline architecture, BFS region growing algorithm, and per-stage parameterization of the generic BFS, see [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md).

---

## General Approach

"Region growing" from seed faces while preserving fit invariants, refitting as necessary as faces are added to the region. Given that vertices lie within epsilon of precisely on the surfaces involved, any seed on a surface should grow to the entire region. The generic BFS algorithm is described in [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md#generic-bfs-hypothesis-generation).

---

## 2.1 Planar fitting math

**Plane fitting** uses **area-weighted normal averaging**: each face's unit normal is weighted by the face's triangle area, then the weighted sum is normalized. The plane distance is the mean of `normal · vertex` over all vertices in the hypothesis. This gives larger faces more influence on the plane orientation, which is more robust than unweighted averaging.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted plane (perpendicular projection).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - (Changed from original plan: the plan described projecting face *edges* onto the hypothesis and checking the enclosed surface region. The centroid-projection approach is simpler and sufficient for validation — it confirms the fitted plane aligns with a real STEP surface at representative interior points. Edge projection would be needed for boundary accuracy validation, which is a stage 3 concern.)

---

## 2.2 Conical and cylindrical fitting math

Cones and cylinders share a single BFS pass. A cylinder is a degenerate cone with half-angle = 0. Both fits are attempted for every seed, and the better one is kept. The seed selection, acceptance test, and BFS parameters are common; the fitting algorithms differ.

### Seed selection

Cone/cylinder seeds use a "fan" topology: initial face $f_i$, neighbor $f_1$, and another face $f_2$ that is a neighbor of $f_i$ but NOT an edge-wise neighbor of $f_1$. This puts $f_i$ in between $f_1$ and $f_2$ for better angular/circumferential coverage than a linear chain. Pairwise qualification checks the normals of $f_1$ and $f_2$ against $f_i$'s normal (not against each other): cross product magnitude must exceed MIN_CROSS_THRESHOLD (0.01) and dihedral angle must be ≤ angular tolerance.

When $f_i + f_1 + f_2$ yield only 5 unique vertices (all triangles sharing edges), this is insufficient for 6-DOF cone fitting. An optional $f_3$ is added — a face adjacent to $f_2$ that doesn't share an edge with $f_i$ or $f_1$ — to bring the count to at least 6 vertices. If no qualifying $f_3$ is found, the 3-face seed is still emitted (sufficient for 4-DOF cylinder fitting).

Seed tolerance multiplier: `seed_vtol_multiplier = 1.0`. Vertices should lie precisely on the conical surface, so the seed vertex tolerance matches the base vertex tolerance.

### Cylinder fitting

**Cylinder fitting** (used for seeding and re-fitting):

The fitting proceeds in three steps:

1. **Axis direction estimation from face normals.** On a cylinder, every face normal is perpendicular to the axis: `n · a = 0`. So the axis direction `a` minimizes $\sum w_i (n_i \cdot a)^2$ subject to $\|a\| = 1$, where $w_i$ is the area of face $i$. This is the eigenvector corresponding to the *smallest* eigenvalue of the 3×3 weighted covariance matrix $M = \sum w_i \, n_i \, n_i^T$. For the three-face seed, this is the full 3×3 eigenvector computation (the cross-product shortcut `a = normalize(n1 × n2)` that works for two-face seeds is not used, since we want the axis estimate that best fits all three faces' normals).

2. **Axis position and radius via 2D circle fitting.** Given the axis direction `a`, choose two orthogonal unit vectors `u`, `w` perpendicular to `a`. Project each vertex onto the 2D plane: `(x_i, y_i) = (v_i · u, v_i · w)`. **Center the 2D coordinates** by subtracting the mean: `x_i' = x_i - mean_x`, `y_i' = y_i - mean_y`. Fit a circle to the centered 2D points using the algebraic least-squares method: fit the linear model $x^2 + y^2 + Dx + Ey + F = 0$ via least squares, then `center' = (-D/2, -E/2)` and `radius = sqrt(center_x'² + center_y'² - F)`. The axis origin in 3D is `(center_x' + mean_x) * u + (center_y' + mean_y) * w` (the component along `a` is arbitrary). Centering is essential for numerical stability: without it, the normal equations matrix becomes ill-conditioned when absolute vertex coordinates are large relative to the arc span (e.g., vertices at (15mm, 12.5mm) spanning 2.5° of a r=2mm cylinder produce sums-of-squares ~1350 that dwarf the r=2 signal, causing wildly wrong radii).
3. **Levenberg-Marquardt refinement.** The initial axis + circle fit provides a good starting point but can have ~0.4° axis error when faces span a narrow arc (<15°). LM refines all 6 parameters simultaneously: `[alpha, beta, qx, qy, qz, radius]` where alpha/beta are tilt angles from the initial axis in two perpendicular directions. **Vertices are recentered** around their centroid before LM to improve numerical conditioning (subtracting the centroid makes the axis origin parameters small, improving Jacobian precision). The residual for each vertex is `||cross(p - q, a)|| - r` (signed distance to cylinder surface). Numerical Jacobian via central differences (eps=1e-8). Uses `levenberg-marquardt` crate v0.14 with `nalgebra` v0.33 types, with `patience=500` for thorough convergence. Falls back to the step-2 estimate if LM produces a negative radius. After optimization, the centroid offset is added back to the refined origin.

The normal covariance matrix `M` can be maintained incrementally during BFS (add $w_i \, n_i \, n_i^T$ when accepting a face), but the eigenvector solve, circle fit, and LM refinement must be recomputed from scratch when re-fitting, since a change in axis direction invalidates the 2D projection. Steps 1-2 are O(n) in vertices with simple arithmetic; step 3 adds ~5-10 LM iterations (typically <1ms per hypothesis).

**Subtleties and edge cases:**
- **Near-parallel seed normals**: Pairwise qualification requires cross product magnitude > `MIN_CROSS_THRESHOLD` (e.g. 0.01, ~0.6° between normals) for each consecutive pair in the seed chain (fi↔n1 and n1↔n2). This ensures each pair contributes meaningful angular information. Note that fi and n2 need NOT have different normals — on a cylinder, faces at the same circumferential position but different axial positions have identical normals. The angular diversity comes from the intermediate face n1 having a different circumferential position.
- **Cones vs cylinders**: A cone has varying radius along the axis. If we fit a cylinder to conical faces, vertex distances near the ends will exceed tolerance, so cone faces will be naturally rejected. Both cone and cylinder fits are attempted for every seed; see the cone fitting subsection below.
- **Large-radius cylinders**: A cylinder with very large radius looks locally planar. Such faces will belong to multi-face planar hypotheses from stage 2.1 (since vertices on a nearly-flat cylinder are nearly coplanar within tolerance). Multi-face planar faces CAN seed cylinders — the angular coverage validation, minimum face count, and other checks prevent bogus large-radius cylinder fits. Stage 2.6 resolves any ambiguity between cylindrical and planar assignments.
- **Cylinder wrapping >180°**: This works naturally with BFS. The convexity check ensures we don't accidentally merge inner and outer surfaces.
- **Multiple cylinders sharing an axis**: E.g. stepped bore holes with different radii. These will form separate hypotheses because vertices at different radii will exceed the vertex tolerance.
- **Seam between cylinder and plane**: Faces at the junction have some vertices on the plane and some on the cylinder. Since the plane vertices aren't on the cylinder surface, the vertex distance check rejects them. This correctly excludes boundary faces.

**Angular coverage validation**:

The problem: any two non-coplanar faces fit *some* cylinder. A bogus seed that doesn't span the cylinder's circumference can grow into a small hypothesis via BFS, consuming faces that should belong to a correct hypothesis from a better seed. The angular coverage check ensures the hypothesis has genuine circumferential support — faces distributed around the cylinder, not clustered on one side.

Algorithm:
1. For each face in the hypothesis, compute the centroid and project it onto the cylinder's angular coordinate θ. Given the cylinder's axis direction `a`, axis origin `o`, and orthogonal basis vectors `u`, `w` perpendicular to `a`: for centroid `c`, compute `radial = c - o - ((c - o) · a) * a`, then `θ = atan2(radial · w, radial · u)`.
2. Sort the θ values (one per face).
3. Compute all N gaps between consecutive sorted θ values, **including the wraparound gap**: `gap_wrap = θ[0] + 2π - θ[N-1]`.
4. The **largest gap** is the "empty arc" — the angular region not covered by any face. The **span** is `2π - largest_gap`.
5. Among the remaining N-1 gaps, the **second-largest gap** must be `≤ span / 3`.

This ensures the hypothesis has at least 3 distinct angular clusters of faces around the cylinder. A bogus cylinder where all faces are side-by-side in a narrow arc will fail: e.g., 5 faces clustered in a 30° arc have span ≈ 30°, but the second-largest internal gap (~7°) may pass. However, such a cluster's centroid validation and vertex tolerance checks should already reject it. The angular coverage check primarily catches the case where a seed grabs faces from 2 sides of a non-cylindrical surface (e.g., two planar faces that happen to fit a cylinder) and BFS adds a few more nearby faces, resulting in a hypothesis with faces clustered in 1–2 angular positions rather than distributed around the circumference.

*Note on applying angular coverage during BFS growth: this was considered but is not useful. Early in BFS, faces are naturally clustered near the seed, so the check would fail prematurely. The vertex distance, convexity, and angular tolerance checks are the correct per-face acceptance criteria during growth. Angular coverage is a global structural validation — it answers "does this completed hypothesis look like a real cylinder?" rather than "should this next face be added?"*

*Note: For partial cylinders (arcs less than 360°), this check still works — it requires 3+ angular clusters within whatever arc is covered, which is the minimum for a reliable cylinder fit. A 90° cylindrical arc with 6 evenly-spaced faces has internal gaps of ~18° and span of ~90°, so second_largest_gap (18°) ≤ span/3 (30°) — passes correctly.*

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the cylinder surface (nearest point: project onto axis, then move radially to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta (the flat triangle's centroid lies inside/outside the curved surface by approximately $h^2 / 8R$ where $h$ is the chord length). For typical tessellations this is small relative to `--surface-tolerance`, but it could become significant for very coarse meshes on tight-radius cylinders.

### Cone fitting

**Cone parameterization:**
A cone with apex $A$, axis direction $\hat{a}$, and half-angle $\theta$ satisfies: for any point $p$ on the surface, the angle between $(p - A)$ and $\hat{a}$ equals $\theta$. Equivalently, in (h,r) coordinates where $h = (p - A) \cdot \hat{a}$ (signed axial distance from apex) and $r = |(p - A) - h\hat{a}|$ (radial distance from axis), the surface satisfies $r = h \tan\theta$.

The signed distance from a point $p$ to the cone surface (for LM residual) is:
$$d = r_i \cos\theta - h_i \sin\theta$$
where $h_i = (p_i - A) \cdot \hat{a}$ and $r_i = |(p_i - A) - h_i \hat{a}|$.

**Two-path fitting algorithm:**

The fitting algorithm detects the cone topology and uses different strategies for pointed vs. truncated cones:

**Pointed cone path** (a vertex shared by ALL seed faces = apex):
1. **Apex** = shared vertex (exact). Validated by checking that non-apex vertex distances from the apex have low coefficient of variation (< 0.5), confirming a conical rather than arbitrary fan pattern.
2. **Axis** = cross product of normal differences. On a cone, $n_i \cdot \hat{a} = \sin\theta$ is constant, so $(n_i - n_j) \cdot \hat{a} = 0$ — all normal differences are perpendicular to the axis. The cross product of two independent normal differences gives the exact axis. Falls back to the largest eigenvector of $\sum \hat{d}_i \hat{d}_i^T$ where $\hat{d}_i = (v_i - \text{apex}) / |v_i - \text{apex}|$ (vertex directions cluster around the axis for narrow cones).
3. **Half-angle** from median of $\arctan(r_i / h_i)$ over all non-apex vertices. The axis is oriented so that $h > 0$ for the majority of vertices (away from apex).
4. **Fixed-apex LM** refinement: 3 parameters (axis tilt $\alpha$, $\beta$ + half-angle $\theta$), keeping the apex pinned. **Vertices are recentered** around their centroid and the apex is shifted accordingly before LM for numerical stability. Uses `patience=500` for thorough convergence. The general 6-parameter LM is underdetermined for pointed cones because the apex vertex has zero residual regardless of parameters, producing near-zero Jacobian rows.

**Truncated cone path** (no shared vertex):
1. **Axis** via multi-strategy estimation:
   - ≥8 faces: smallest eigenvector of area-weighted normal covariance $M = \sum w_i n_i n_i^T$ (same as cylinder).
   - <8 faces with a vertex shared by most faces: direction from that vertex to the centroid of remaining vertices.
   - <8 faces, no dominant vertex: cross product of normal differences (as in pointed path).
2. **Apex and half-angle** from linear regression of vertex projections in (h,r) space: $r_i = m \cdot h_i + b$ where $h_i = (v_i - \text{centroid}) \cdot \hat{a}$ and $r_i$ is radial distance. Then $\theta = \arctan|m|$ and apex offset along axis $= -b/m$.
3. **Full LM** refinement: 6 parameters ($\alpha$, $\beta$, apex $x,y,z$, half-angle). Analytic Jacobian. **Vertices are recentered** around their centroid before LM (apex parameters are also shifted), then the centroid offset is added back to the refined apex. Uses `patience=500` for thorough convergence.

**Fit selection (cone vs. cylinder):**

Both cone and cylinder fits are attempted for every seed. Max vertex error and max surface (centroid) error are computed for each; the fit with lower max vertex error is selected. The rejection criteria (half-angle bounds, apex distance check, angular coverage) are applied independently.

### Acceptance test (cones and cylinders)

1. **Minimum face count**: cones require MIN_CONE_FACES, cylinders require MIN_CYLINDER_FACES.
2. **Half-angle bounds**: reject cones with $\theta > 80°$ (nearly planar).
3. **Angular coverage**: faces must be distributed around the cone/cylinder circumference (algorithm described above).

*Note: A normal-axis consistency check was previously used (rejecting hypotheses where the variance of $n_i \cdot \hat{a}$ exceeded 0.0005), but was removed. It's acceptable for cone hypotheses to be detected on spherical patches — the surface selection phase (stage 2.6) picks the sphere hypothesis over conical hypotheses because spheres cover more area.*

### BFS parameters

Cone BFS uses `angular_tol * 2.0` because tessellated cones with $N$ circumferential divisions have inter-strip dihedral angles $\approx 360°/N$. For $N=20$ (typical), that's ~18°, which exceeds the default angular tolerance of 17.5°. Vertex-to-surface distance is the primary discriminator for cones, not angular similarity.

### Performance notes

1. **Filter by cylindrical hypothesis**: Skip faces already committed to cylindrical hypotheses when building apex candidates, collecting seed faces, and during BFS expansion.
2. **Early half-angle rejection**: After initial cone fit but before seed tolerance check and BFS, reject cones with half-angle outside [2°, 85°].
3. **Do NOT filter by spherical hypothesis**: Cone faces can have spurious spherical hypotheses (e.g., `block_with_conical_hole`: all 94 cone faces get spherical hypotheses). Only cylindrical filtering is safe because cones and cylinders are geometrically distinct, while cones and spheres overlap locally.

- With the `--compare` flag:
    - Same centroid-projection validation as cylinder/sphere.

---

## 2.3 Spherical fitting math

**Sphere fitting** (used for seeding and re-fitting):
- Given vertices on a sphere: $|v - c|^2 = r^2$, expand to $v \cdot v - 2 v \cdot c + |c|^2 = r^2$.
- Rearrange: $v \cdot v = 2 v \cdot c - |c|^2 + r^2 = 2 v \cdot c + k$ where $k = r^2 - |c|^2$.
- This is linear in unknowns $(c_x, c_y, c_z, k)$. Solve via least squares: $A x = b$ where $A_{i} = [2v_{ix}, 2v_{iy}, 2v_{iz}, 1]$ and $b_i = v_i \cdot v_i$.
- After solving: center $= (c_x, c_y, c_z)$, radius $= \sqrt{k + |c|^2}$.
- **Gauss-Newton refinement** (when ≥6 vertices): After the algebraic fit, iteratively minimize the geometric distance $f_i = \|p_i - c\| - r$. The Jacobian is $J_i = [-(p_i - c)/\|p_i - c\|, -1]$. Solve the normal equations $J^T J \delta = -J^T f$ at each iteration (4×4 system via Gaussian elimination). Up to 20 iterations, converges when $\|\delta\| < 10^{-14}$. This refines the algebraic fit to minimize the true vertex-to-surface distance rather than the algebraic residual.

**Maximum radius**: `max_sphere_radius`. With solid-angle coverage validation and surface-tolerance validation during BFS, the radius limit no longer needs to be tight — those checks prevent the pathological growth patterns that the original 10× bounding box diagonal limit was designed to catch. A large value like `bounding_box_diagonal * 1000` suffices as a numerical guard rail, preventing sphere fits with absurd radii (where floating-point precision of vertex coordinates could make flat faces appear spherical) while not rejecting any plausible real-world spherical feature. Compute the bounding box once at stage 2 entry.

**Solid-angle coverage validation** (applied to the completed hypothesis, analogous to angular coverage for cylinders):

The problem: sphere BFS can grow along cylinder fillet surfaces, because locally adjacent cylinder faces fit on a sphere within vertex_tolerance. This produces oversized sphere hypotheses spanning a thin strip around a cylinder fillet (observed on `onshape_rounded_cube`). The analog of the cylindrical angular coverage test is a **solid-angle coverage test**: faces on a genuine spherical patch produce centroid-to-center directions that span 3D, while faces on a cylinder fillet produce directions that span only a 1D strip.

Algorithm:
1. For each face in the hypothesis, compute the unit vector from sphere center to face centroid: $d_i = \text{normalize}(\text{centroid}_i - \text{center})$.
2. Compute the area-weighted 3×3 covariance matrix of these unit vectors: $C = \sum w_i \, d_i \, d_i^T$ where $w_i$ is the face's mesh area.
3. Compute eigenvalues of $C$. Let $\lambda_1 \geq \lambda_2 \geq \lambda_3$.
4. Require $\lambda_3 / \lambda_1$ to exceed a minimum threshold (e.g., 0.01). If it fails, the hypothesis directions are nearly coplanar (a strip on the sphere, characteristic of cylinder fillet growth) or nearly collinear (a narrow band), and the hypothesis should be discarded.

*Note: Fillet faces lie along a 1D arc on the sphere's direction space, producing $\lambda_3 \approx 0$. Genuine spherical patches (hemispheres, domes, spherical pockets) have centroid directions spanning 3D, so all eigenvalues are substantial. For partial spherical caps, the directions span at least 2D; very small caps may fail but such small patches are unlikely to be meaningful spherical features.*

*Note: Like cylindrical angular coverage, this should NOT be applied during BFS growth — early in BFS, faces are clustered near the seed and naturally appear low-dimensional. It is a global structural validation for completed hypotheses.*

**Distinguishing spheres from cylinders:**
- Faces with cylindrical hypotheses are NOT excluded — they may legitimately belong to both a cylinder and a sphere (e.g., equator band). Stage 2.6 disambiguates.
- On cylindrical surfaces, sphere fits from vertex neighborhoods with nearly-parallel normals tend to produce unreliable results (very large or negative $r^2$), which are caught by `max_sphere_radius` and the $r^2 > 0$ check in sphere fitting. The solid-angle coverage test provides additional protection: faces along a cylinder strip span only 1D in direction space, failing the eigenvalue ratio test.
- Vertex-based seeding with the minimum-3-face requirement and angular diversity naturally avoids seeds from planar regions while allowing seeds from curved regions.

**Preventing false positives on chamfers/bevels:**
- On a chamfered or beveled cube, corner and edge chamfer faces can seed sphere fits. However, these faces are flat — the sphere fit produces a degenerate radius vastly larger than the model. The `max_sphere_radius` constraint catches this even at 1000× bounding box diagonal: for a 14mm model with vertex_tolerance=1e-5, a degenerate sphere fit would have R ≈ 112 km, far exceeding `14mm × 1000 = 14m`.
- The centroid validation check also catches spurious fits where the algebraic fit succeeds but the sphere doesn't actually pass through the face centroids.
- The solid-angle coverage test may also reject chamfer fits where few faces span a limited set of directions.

**Test models:**
- Simple sphere (ccad): full convex sphere, radius 10 → 1 convex spherical hypothesis.
- Hemisphere (ccad): top half of sphere + flat base → 1 convex spherical hypothesis + 1 multi-face planar.
- Spherical pocket (ccad): block with concave hemispherical cavity → 6 planar + 1 concave spherical.
- Ball on cylinder (ccad): sphere atop cylinder stalk → 1 convex spherical + 1 convex cylindrical + 1 planar.
- Sphere (onshape): full sphere, fine tessellation.
- Dome/hemisphere (onshape): hemisphere with flat base.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted sphere surface (nearest point: normalize the centroid-to-center vector and scale to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: same centroid-to-surface sagitta caveat as for cylindrical hypotheses (stage 2.2). For spheres the sagitta is approximately $h^2 / 8R$ where $h$ is the chord length of the facet.

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. Also consider toroidal surfaces (donuts, fillets)—they're common in CAD and combine characteristics of cylinders and spheres. Note: torus faces locally fit spheres, producing sphere hypotheses on torus surfaces (e.g., pipe_elbow). Stage 2.6 will need torus hypothesis support to resolve this.*

---

## 2.4 Toroidal fitting math

**Torus parameterization:**
A torus with center $C$, axis $\hat{a}$, major radius $R$, minor radius $r$ satisfies: for any point $p$, let $v = p - C$, $\text{axial} = v \cdot \hat{a}$, $\text{radial} = |v - \text{axial} \cdot \hat{a}|$. Then:
$$d = \sqrt{(\text{radial} - R)^2 + \text{axial}^2} - r$$
is the signed distance from $p$ to the torus surface (0 on surface, positive outside tube, negative inside).

**Fitting via medial axis / tube center method:**

Direct nonlinear fitting of 7 torus parameters ($C$, $\hat{a}$, $R$, $r$) is numerically fragile because of the high DOF and the coupling between parameters. The **medial axis approach** linearizes the problem by exploiting the geometric relationship between surface normals and the tube center.

Key insight: for a surface with constant positive principal curvature $\kappa_1 = 1/r$ in one direction, the **tube center** (center of the osculating circle) at each point is $k_i = p_i + r \cdot n_i$ where $n_i$ is the outward surface normal. For a torus, all tube centers lie on the **major circle** — a circle of radius $R$ centered at $C$ in the plane perpendicular to $\hat{a}$. This reduces torus fitting to: (1) estimate $r$, (2) compute tube centers, (3) fit a 3D circle.

**Step 1: Estimate minor radius $r$ from normal-line intersections.**
- For pairs of adjacent (non-parallel) faces in the candidate region, compute the closest approach of their normal lines: $L_1: p_1 + t \cdot n_1$, $L_2: p_2 + s \cdot n_2$.
- Reject pairs with nearly parallel normals ($|n_1 \cdot n_2| > 0.95$).
- The closest-approach distance approximates 0 (for a torus, normal lines converge at the tube center), and the intersection parameter gives $r \approx t \approx s$ (for outward normals on a convex torus).
- Collect many estimates, take the **25th percentile** for robustness. On a torus, normal-line intersections measure a mix of tube curvature ($\approx r$, the minor radius) and azimuthal curvature ($\approx R \pm r$). The tube-curvature estimates are the smallest values, so a low percentile avoids contamination from azimuthal estimates. (Median would pick $\approx R - r$ for concave patches.)

**Step 2: Compute tube center points.**
- For each vertex $p_i$ with outward normal $n_i$: $k_i = p_i + r \cdot n_i$ (convex torus) or $k_i = p_i - r \cdot n_i$ (concave torus).
- These points should cluster on the major circle.

**Step 3: Classify spine geometry (detect degeneracies).**
- Compute PCA of the $k_i$ points. Eigenvalue pattern:
  - $\lambda_1 \approx \lambda_2 \gg \lambda_3$: points form a circle → **torus** (proceed).
  - $\lambda_1 \gg \lambda_2 \approx \lambda_3$: points form a line → actually a **cylinder** (delegate to cylinder detection).
  - $\lambda_1 \approx \lambda_2 \approx \lambda_3 \approx 0$: points cluster at a single point → actually a **sphere** (delegate to sphere detection).

**Step 4: Fit 3D circle to tube center points.**
- PCA of $k_i$: smallest eigenvector = circle plane normal = **torus axis direction** $\hat{a}$.
- Project $k_i$ into the circle plane using an orthonormal basis $(u, v)$ perpendicular to $\hat{a}$.
- 2D algebraic circle fit in the projected plane → center $(u_c, v_c)$ and radius $R$.
- Convert back to 3D: torus center $C = \text{mean}(k_i) + u_c \cdot u + v_c \cdot v$, major radius $R$.

**Step 5: LM refinement.**
- Optimize $[C_x, C_y, C_z, \alpha, \beta, R, r]$ (7 parameters: center 3, axis tilt 2, major radius 1, minor radius 1).
- Residual: $d_i = \sqrt{(\text{radial}_i - R)^2 + \text{axial}_i^2} - r$.
- LM patience 500, minimum 7 data points. Analytic Jacobian with closed-form partial derivatives for all 7 parameters (center, axis tilt, major and minor radii).

**Alternative fitting: SOR axis method (`fit_torus_sor`).**
When the tube center method fails (e.g., too few tube centers above threshold or quality gate failure), a surface-of-revolution axis estimation is used as fallback. Estimates torus axis from PCA of vertex positions (smallest eigenvector of position covariance matrix) — for a torus, axial variance is ~$r^2/2$ while perpendicular variance is ~$(R^2 + r^2/2)/2$, so the axis has minimum positional variance for any non-self-intersecting torus ($R > r$). Note: face-normal covariance does NOT give the axis as its smallest eigenvector for a torus — unlike cylinders where normals are perpendicular to the axis, torus normals have significant axial components. Center found via Kasa circle fit on 2D projections of vertices perpendicular to axis. Projects vertices to (h,r) coordinates, fits a shifted circle to deduce R and r. Same LM refinement follows. `fit_torus` tries **both** methods (tube center and SOR) and keeps whichever achieves the lowest maximum vertex-to-torus distance.

**Multi-axis LM rescue:**
If the initial geometric fit (tube center or SOR) does not meet vertex_tolerance, a multi-axis Levenberg–Marquardt rescue is attempted. Five axis candidates are tried: the geometric estimate, the smallest eigenvector of the face-normal covariance matrix, and the three coordinate axes [1,0,0], [0,1,0], [0,0,1]. For each axis candidate:
1. Project vertices to the plane perpendicular to the axis and fit a 2D algebraic circle (Kasa method, 3×3 linear system) to estimate center position and major radius $R$.
2. Estimate minor radius $r$ from the residuals of the 2D circle fit.
3. Run LM refinement from this initial estimate.
4. Track the best result across all axis candidates. Accept if the best achieves vertex_tolerance.

**Deterministic iteration.**
All vertex set iterations (HashSet) are sorted before use to ensure deterministic results regardless of hash ordering. Without this, the fitting results can vary between runs due to floating-point non-associativity.

*Comment: The medial axis / tube center method is far more robust for fillet-sized torus patches than direct 7-parameter nonlinear fitting. By reducing the problem to (1) estimate one scalar r, (2) fit a 3D circle to the tube centers, it separates the easy part (circle fitting is linear algebra) from the hard part (estimating the minor radius). The spine classification step (Step 3) also provides a natural way to detect degenerate cases where the "torus" is actually a cylinder or sphere — this happens when the fillet radius approaches the edge length or when the surface is nearly spherical.*

---

## Algorithm Reference

This section documents key algorithmic techniques drawn from computational geometry research that inform the surface fitting approach. These are referenced by the stage descriptions above.

### Unified Surface-of-Revolution (SOR) Framework

All surfaces of revolution — cylinder, cone, sphere, torus — share a common structure: they are generated by rotating a 2D profile curve around an axis. This means they can all be fitted using a unified approach:

1. **Estimate the axis of revolution** from the surface normals: for any SOR, surface normals are perpendicular to the axis (cylinder), make a constant angle with the axis (cone), point toward the axis (sphere pole behavior), or some combination. The axis can be estimated as the smallest eigenvector of the area-weighted normal covariance matrix $M = \sum w_i n_i n_i^T$.

2. **Reduce to (h,r) coordinates**: for each vertex, compute $h = (p - \text{centroid}) \cdot \hat{a}$ (axial position) and $r = |(p - \text{centroid}) - h\hat{a}|$ (radial distance from axis).

3. **Fit candidate profile curves** in the (h,r) plane:
   - **Cylinder**: $r = R$ (constant) — horizontal line. 1 parameter: $R$.
   - **Cone**: $r = m \cdot h + b$ (linear) — sloped line. 2 parameters: $m, b$.
   - **Sphere**: $h^2 + (r - 0)^2 = R^2$ — circle centered on the axis. 2 parameters: $h_0, R$ (center offset and radius). But note: sphere center is ON the axis, so this is $h^2 + r^2 = R^2$ in centered coordinates.
   - **Torus**: $(r - R)^2 + h^2 = r_{\text{minor}}^2$ — circle offset from axis. 2 parameters: $R, r_{\text{minor}}$.

4. **Pick the best profile** by RMS residual. Use BIC or AIC to penalize model complexity (cylinder has 1 DOF, cone has 2, torus has 2).

This framework naturally handles ambiguous cases: a very gentle cone looks like a cylinder (slope $\approx 0$); a sphere with very large radius looks like a plane ($R \to \infty$); a torus with large major radius looks like a cylinder ($R \gg r$).

The current implementation uses separate fitting algorithms for each surface type (stages 2.2–2.4), which is appropriate for incremental development. A future refactoring pass (see Implementation Phases in DEVELOPMENT_PLAN.md) will unify these into the SOR framework to improve accuracy on ambiguous cases and simplify the codebase.

### Gaussian Map (Normal Covariance) Classification

The distribution of surface normals (the "Gaussian map") directly reveals the surface type. For a patch of mesh faces, compute the area-weighted covariance matrix of the normals:
$$C = \sum_i w_i (n_i - \bar{n})(n_i - \bar{n})^T$$
where $w_i$ is face area and $\bar{n}$ is the weighted mean normal. The eigenvalues $\lambda_1 \geq \lambda_2 \geq \lambda_3$ classify the surface:

- **Plane**: $\lambda_1 \approx \lambda_2 \approx \lambda_3 \approx 0$ — all normals identical.
- **Cylinder/Cone**: $\lambda_1 \approx \lambda_2 \gg \lambda_3$ — normals lie on a great circle (cylinder) or small circle (cone).
- **Sphere**: $\lambda_1 \approx \lambda_2 \approx \lambda_3$ (all nonzero) — normals distributed isotropically.
- **Torus**: $\lambda_1 > \lambda_2 \gg \lambda_3$ — normals form a band pattern.

The cylinder axis estimation already uses this: the axis is the eigenvector of the uncentered normal covariance $M = \sum w_i n_i n_i^T$ corresponding to the smallest eigenvalue (the direction perpendicular to all normals). Extending this to pre-classify surface type before fitting would improve seeding efficiency.

### Medial Axis / Tube Center Method

For any surface with a constant principal curvature $\kappa_1 = 1/r$ in one direction, the **medial axis** (locus of centers of curvature) can be computed directly from the surface points and normals:
$$k_i = p_i + r \cdot n_i$$

The geometry of the medial axis reveals the surface type:
- **Cylinder**: medial axis is a line (the cylinder axis) — $k_i$ points are collinear.
- **Sphere**: medial axis is a point (the center) — $k_i$ points cluster.
- **Torus**: medial axis is a circle (the major circle) — $k_i$ points lie on a circle.
- **Cone**: medial axis is a point (the apex) for one direction of curvature.

This is used in stage 2.4 (torus fitting): estimate $r$ from normal-line intersections between adjacent faces, compute tube centers $k_i$, fit a 3D circle to get the major circle, and extract the torus parameters.

### Quadric Fitting

Any quadric surface (plane, sphere, cylinder, cone, ellipsoid, hyperboloid, paraboloid) can be written as:
$$Q(x, y, z) = Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0$$

Fitting: construct the $N \times 10$ feature matrix and find its null space via SVD (smallest singular vector). The matrix form $x^T \mathbf{A} x + b^T x + c = 0$ encodes the surface type in the eigenvalues of $\mathbf{A}$:
- **Sphere**: $\mathbf{A} = \lambda I$ (three equal eigenvalues).
- **Cylinder**: one eigenvalue is 0 (axis direction), other two equal.
- **Cone**: three nonzero eigenvalues with $\det(\mathbf{A}) \approx 0$ (passes through apex).
- **Plane**: $\mathbf{A} = 0$ (all eigenvalues zero).

Quadric fitting is fast (SVD of an $N \times 10$ matrix, $O(N)$) and provides a good initialization for specialized fitters. It could be used as an alternative initialization for cone fitting (stage 2.2).
