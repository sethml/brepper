# Surface Fitting Algorithms

This document describes the surface fitting math and validation algorithms used in brepper's Stage 2 pipeline. Each section covers the fitting math and validation for one surface type.

> **Orientation (post stage-2 rewrite).** The per-type **fitting math and
> validation** below (the `fit_surface`/`vertex_distance`/convexity content)
> is current — those fitters now implement the trimmed `surface::SurfaceType`
> trait and are driven generically by the whole-region DOF ladder
> (`classify.rs`) and scoped competition (`compete.rs`). The **BFS region
> growing / seeding / acceptance-test / seeding-mask** orchestration this
> document still refers to in places was the old pipeline and **no longer
> exists** (removed in Rewrite Roadmap step 7); read those mentions as
> historical. See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) §2 for the
> current partition-then-classify architecture.
>
> **Closed-form initialization (Phase 1.1/1.2, 2026-07-31).** A newer module,
> `src/stage2/initialize.rs`, fits each analytic type from a handful of facets
> in closed form — see "Closed-form initialization from facet normals" below.
> It is **not wired into the pipeline**: the earn loop that was to consume it
> (Phase 1.2) was measured and not landed, so today the module is exercised only
> by its unit tests. Its math and its measured failure modes are recorded here
> because they are properties of the geometry, not of that attempt.

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

## Precision and conditioning (measured 2026-07-28)

What limits a fit is not the fitter but the data and the support. Three facts,
all measured on the corpus, that the fitting math has to respect.

**The accuracy floor is a float32 quantum, and it scales with the model.** The
corpus ships *binary* STLs, so vertices are float32 and the floor on vertex
accuracy is `max|coordinate| × 1.19e-7`: 2.76e-6 mm for rounded_cube (max coord
23.1, measured stage-1 vertex distance 1.9e-6) and **1.96e-5 mm for
stepped_block_toroidal** (max coord 164, measured 1.23e-5). So "a mesh vertex
lies on the true CAD surface to ~1e-6 mm" is *not* a constant — it is twenty
times worse on a 164 mm part than on a 23 mm one, and a part exported 200 mm
off-origin carries ~10× worse data than one centred. Nothing downstream recovers
it. Expressed in quanta, committed analytic fits land at 0.3–0.5 (median) across
a 7× range of quanta, which is a scale-free way to say "as well as the data
allows".

**A short arc cannot determine a free surface of revolution, by any method.**
For an arc of half-angle θ on radius r the sagitta is ≈ rθ²/2. Over three
2.5° strips (7.5° of arc) on r = 0.5 that is ~0.001 mm against a 2e-5 mm
quantum, so the radius comes out ~2% wrong — 0.01 mm, or 200× vertex tolerance.
No fitting method repairs this; only more arc does. This is why an initial fit
taken from a small patch cannot be trusted, and why growth criteria must not
depend on the fit already being accurate.

**Constraining the axis collapses the problem.** Measured on the three fillet
cylinders `stepped_block_toroidal` fine fails to find: fixing *only the axis
direction* and fitting a 2D circle in the perpendicular plane recovers
**r = 0.99998, 0.99999 and 0.50008** (true 1.0, 1.0, 0.5) from arcs of 65°, 60°
and 40°, with maximum residuals of 1.0e-5, 1.0e-5 and 5.0e-6 mm — five to ten
times *under* vertex tolerance. The free 5-DOF fit fails on the **axis**, not the
radius. An adjacent committed surface supplies exactly that constraint (see
`locus.rs` and DEVELOPMENT_PLAN, "Rewriting stage 2.3(b)").

**Corollary for identifiability.** Residual cannot police over-parameterised
fits: the three spurious tori on that model fit their 7–8 faces to 0.52/0.72/0.87
quanta against genuine surfaces at mean 0.39, max 0.56 — overlapping, because a
7-DOF torus on 7 faces is essentially exactly determined. Guards must ask whether
the support *pins* the parameters (angular sweep, conditioning), never how small
the residual is.

---

## Closed-form initialization from facet normals (measured 2026-07-31)

`src/stage2/initialize.rs`. Each analytic type fitted from a handful of facets
with no iteration, exploiting the fact that **facet normals carry exact
information** where a free least-squares fit on tiny support does not (see
"A short arc cannot determine a free surface of revolution" above). Unit-tested
to machine precision on synthetic strips. Not currently called by the pipeline.

### The constructions

- **Cylinder** — axis $= n_1 \times n_2$, exact from two facets and needing no
  angular spread at all; then radius and position from the existing
  `fit_circle_for_axis` over the patch's vertices.
- **Sphere** — the line through a facet's **circumcenter** along its normal
  passes exactly through the centre (the circumcenter is equidistant from three
  points of the sphere, so it lies on the axis of their circle). Two facets give
  the centre as the common-perpendicular midpoint; the radius is the mean vertex
  distance. The *centroid* does not have this property — only the circumcenter.
- **Cone** — see below; this one is not a single construction.
- **Torus** — the "bilinear" meridian-coplanarity system is **linear in Plücker
  coordinates**. With the axis line's moment $m = a \times d$ and
  $q_i = p_i \times n_i$, the constraint $(p_i - a)\cdot(n_i \times d) = 0$
  becomes $q_i \cdot d + n_i \cdot m = 0$ — one linear equation per facet in six
  unknowns, solved as a null space with no alternation and no starting guess.
  Geometrically it is just "every normal line of a surface of revolution meets
  the axis". The system must be normalized to the patch centroid and extent
  first: the $q$ block carries length units and the $n$ block does not, so on a
  patch far from the origin the null space drowns in the scale mismatch. Then
  $(R, r)$ from the existing `fit_circle_2d` on the $(\rho, z)$ profile.

### Why the cone needs a different parameterization

$(\text{apex}, \text{axis}, \text{half-angle})$ puts a point that runs off to
**infinity** as a patch flattens toward a cylinder at the centre of the estimate,
and any estimator for an unbounded quantity explodes near its degenerate case.
Measured on three adjacent facets of an 18.4° cone it returned a half-angle of
**87°** with the apex halfway up the model, and then validated **none** of the
three facets it had been built from.

Parameterize by the **taper** $dρ/dz = \tan\alpha$ instead — the cone's analogue
of curvature. `taper = 0` is a *cylinder*: a finite value in the middle of the
range rather than `apex = ∞` at its edge, exactly the reason circle fitters use
$\kappa = 1/r$ rather than $r$. Both shape parameters then come from bounded
local quantities — axis from the centred normal covariance, $\sin\alpha =
\overline{|n\cdot d|}$ — and only the *position* is solved, as a linear least
squares in $(u_1, u_2, z_A, m)$ from
$|r|^2 - \tan^2\!\alpha\,z^2 = 2r\cdot u - 2\tan^2\!\alpha\,z\,z_A - m$.

Two routes that do **not** work for a cone, both measured:

- **Normal-line concurrence** (the torus construction) is rank-deficient on a
  cone, and the degenerate case is what a real cone tessellation *is*: a cone is
  ruled, so one row of quads spans its whole height, every facet centroid sits at
  the same axial position, and their normal lines then all meet the axis at the
  **same point**. Any line through that point solves the system. On a torus the
  meeting point slides with the tube angle, which is why the identical solve is
  well posed there. "Surface of revolution" is not the enabling property —
  *variation in where the normal lines cross the axis* is.
- **Uncentred normal covariance** flips above 35.26° — see the Gaussian-map
  section.

### The limit of all of it: normals are chord planes

A facet's plane is a **chord** plane, not the tangent plane. For a strip quad
spanning two rulings the lean works out to
$n\cdot d = \sin\alpha\,\cos(\Delta\theta/2)/\sqrt{1 - \sin^2\!\alpha\sin^2(\Delta\theta/2)}$,
an $O(\Delta\theta^2/8)$ bias — predicted 1.5e-3 relative at
$\Delta\theta = 0.11$ rad and measured as 9.985° against a true 10.000°.

That bias is small against a *wrong answer* but large against a *tolerance*: at a
15.65° tessellation step it is 0.93%, roughly **200× the vertex tolerance**
fits are validated at. So a closed form built on facet normals cannot reach the
validation bar on a real mesh. Only quantities read from **vertices** — which lie
on the surface exactly — can. The practical rule:

> Take directions from normals, magnitudes from vertices, and treat a
> normal-derived shape parameter as a **seed for least squares**, never as the
> answer.

### Ruling alignment: when a facet normal is exact at all

A facet normal on a ruled surface is perpendicular to the axis **iff two of its
vertices share a ruling**. Project the vertices along the axis: they land on a
circle, the facet's plane is parallel to the axis iff those projections are
collinear, and three points on a circle are collinear only when two coincide.
So a triangle straddling three rulings tilts its normal off perpendicular
however perfect the underlying surface is.

Measured over the corpus's 12 cylinder-bearing models, 12 275 cylinder facets:
facets with an edge on a ruling are perpendicular to a worst case of
**0.0006°**, while the **25.6%** that straddle reach **1.0528°** — and every
non-zero perpendicularity in the corpus comes from a straddling facet. Quad
fusion is a 23× discriminator but not a guarantee (2.7% of fused quads straddle
against 62.3% of triangles).

**Tapering surfaces make alignment impossible in places.** A mesher holding
chord length roughly constant drops circumferential divisions as the radius
shrinks: `cone_cylinder` runs 32 divisions at its base ring, 23 through the
middle, 29 where the cylinder starts. Where consecutive rings disagree, *no*
facet can be ruling-aligned — no two vertices share a $\theta$ — and the strip
cannot fuse into quads. Where they agree, the facets fuse and their normals
agree to 7 digits.

### Conditioning: the amplification law

The axis error from a two-facet construction is
**$\text{perpendicularity error} / |n_1 \times n_2|$**, verified to 0.1% on
`pipe_elbow_10` Parasolid medium (a perpendicularity max of 0.0280° over a
conditioning p50 of 0.0271 predicts 1.033°; measured 1.0342°, constant across
all 348 pairs).

Adjacent facets differ in normal by the exporter's angular step, so
$|n_1 \times n_2| \approx \sin(\text{step})$ — p50 0.0432 on a fine mesh, which
is $\sin(2.48°)$ and independently recovers Onshape's 2.5° setting. **The
amplification is therefore ~23× on a fine mesh and gets worse as the mesh
refines**, which inverts the natural assumption that denser meshes fit better.
Pair the angularly *most distant* facets, never neighbours. A 90° reading is the
fully degenerate pair of two facets on one ruling strip: no information, rather
than a wrong answer.

**Selection is not collection.** Spread is the right rule for choosing which
facets to *initialize from*; it is the wrong rule for choosing which faces to
*collect*, because "most different normal" also describes a face on the
neighbouring surface. Collecting by spread split every fillet of
`rounded_cube_offset_occt_medium` (12 cylinders → 24) and invented two spheres on
`cone_tangent_planes_occt_medium`. Curvature distinguishes surfaces; normal
direction does not — but curvature is not sufficient either, since a rounded
cube's fillet radius equals its corner-sphere radius. See DEVELOPMENT_PLAN,
"Phase 1.2 as attempted".

---

## Cone initialization: the options not yet built (surveyed 2026-08-06)

The cone is the worst-initialized surface in the pipeline, and unlike the torus
its problem is not a single fixable defect. This section records what fails, what
the literature says is and is not possible, and the three routes that remain
open — none of them implemented. Nothing here is a recommendation; it is the
option set, with each option's point/facet cost and its known weakness.

### The measured starting point

On `simple_cone` (true half-angle 9.4623°, apex at (0,0,60)), a two-quad seed
whose six vertices lie in two rings of three spanning 31.3° of arc:

| quantity | value |
|---|---|
| reference parameters' residual on those six vertices | 4.1e-7 mm = 0.01× vtol = 0.40 quanta |
| `fit_cone`'s half-angle | **2.54°** |
| `axis_normals_consistent` on that fit | correctly refuses it |
| the fallback cylinder's residual on its own seed | 1.1e-2 mm = **222.7× vtol** |

So the information is present and exact, and every existing route misses it. The
routes and why they decline, at this support size:

| route | requirement | on a 2-facet seed |
|---|---|---|
| `cone_axis_from_centred_normals` (`initialize.rs`) | centred normals spanning two directions | rank 1 → declines |
| `estimate_axis_from_normal_diffs` (`conical.rs`) | ≥3 normals | declines |
| `compute_truncated_cone_axis`, <8 faces | — | apex→centroid *guess* |
| `fit_cone_normal_line_rescue` | normal-line pairs | narrow arcs degenerate |

Fraction of cone facets that have *some* conforming seed set which grows to cover
their surface, measured with `--viz-seeds`: `simple_cone` 4/12, `cone_cylinder`
7/12, `nosecone` 10/12.

### What the literature rules out

Andrews & Séquin (CAD&A 2013) state that no linear least-squares method is known
for fitting circular cones or cylinders from **point positions alone**, citing
Petitjean's survey; linear methods exist only for point clouds *with normals*
(Lukács, Martin & Marshall, ECCV 1998). The obstruction is structural rather than
a gap in the literature: a right circular cone's quadric form must be
`M = d dᵀ − cos²α·I`, i.e. `M` must have a **repeated eigenvalue**, and that
condition is neither linear nor quadratic in the quadric coefficients, so it
cannot be imposed on a linear fit.

Read that as a bound on what to look for. "Analytic, points-only, six points,
no structural assumption" is not available.

### Option A — kinematic scaling field for the apex

Fit the linear velocity field `v(p) = γp + a` and require it tangent to the
surface, `v(p)·n(p) = 0`. Each oriented facet gives one linear equation

```text
γ (p·n) + a·n = 0
```

in four unknowns, homogeneous — a null-space solve from **3 facets**. The apex is
`A = −a/γ`. Then translate the apex to the origin and fit the centred quadric
`Σ cᵢⱼ pᵢpⱼ = 0` (Taubin) from **5 points**; `M`'s eigenvalues are
`(sin²α, −cos²α, −cos²α)`, so the odd eigenvector is the axis and
`tan²α = λ_single / −λ_repeated`. Closed form throughout, no iteration.

**Why this is the interesting one.** Divide the equation by γ and it reads
`n·(p − A) = 0` — *exactly* Phase 1.1's `nᵢ·A = nᵢ·pᵢ`, the relation already in
`initialize.rs`. What differs is the parameterization, and that is precisely the
objection Phase 1.1 recorded when it rejected the apex solve: "the apex runs off
to infinity as a patch flattens toward a cylinder, and an estimator for an
unbounded quantity explodes near its degenerate case". In `(γ, a)` coordinates
**γ → 0 is the cylinder** — a regular interior point of the parameter space, not
an edge — and the paper's own treatment is to read γ below a floor as "fit a
cylinder instead". This is the same Pratt/Taubin move DEVELOPMENT_PLAN already
adopts for the bent tube (`k = 1/R`, "no infinity anywhere"), applied to the apex
rather than the spine curvature.

**Cost and weakness.** Needs 3 facets, so a two-quad seed is out — with two
normals no normals-based method can work, and that is structural. The plane
relation is exact for a facet spanning two rulings and carries the
`O(Δθ²/8)` chord-plane lean otherwise (measured 1.5e-3 relative at Δθ = 0.11
rad), so it is near-exact rather than exact on a straddling tessellation. Ample
for an initializer.

### Option B — general quadric, then eigen-projection

Fit an unconstrained quadric by linear least squares (10 coefficients, 9 DOF, so
**≥9 points**), then project the 3×3 form onto `d dᵀ − cos²α·I` by averaging its
two closest eigenvalues. Points only, no normals, no structural assumption — the
closest thing to the stated preference that exists.

**Cost and weakness.** Algebraic rather than geometric distance, so it is biased
and is an initializer only. Nine points is reachable at modest support (a 4-quad
vertex star has 9 vertices, a 4-quad chain 10) but rules out the 6-vertex chain.
The existing `### Quadric Fitting` note below already flags this as a possible
cone initialization; this is what it would cost.

### Option C — two circles from the ring structure

Two ruling-aligned quads give six vertices in two rings of three, one per rim.
Three points in space determine a unique circle, so each rim yields a centre, a
radius and a normal; the line through the two centres is the axis and the two
radii at two heights give the half-angle by similar triangles. Exact, vertices
only, non-iterative, and it works from **2 facets**.

Measured on the seed above: axis **0.0002°** off, half-angle **0.0001°** off,
the two radii recovered to 0.1–0.2× vtol, against a mesh error of 4.1e-7 mm.
Conditioning is comfortable rather than marginal — the rings span 31.3°, and the
sensitivity `2ε/α²` gives a radius uncertainty of 1.1e-5 mm, under a quarter of
vtol.

**Cost and weakness.** Needs the ruling structure the other options do not: two
rings with ≥3 points each, i.e. ruling-aligned facets. Phase 1.0 measured 25.6%
of corpus cylinder facets straddling three or more rulings, so this cannot be the
only route. Ring identification needs no provisional axis — the shared edge
between two quads in a strip *is* a ruling, so its endpoints are one per ring and
walking the strip labels the rest combinatorially.

### Considered and rejected

- **Solving the 6-DOF system from exactly six points.** Exactly determined, so
  finitely many right circular cones pass through them — but a polynomial system
  with multiple roots and no least-squares character. It is also the configuration
  where a *second exact cone* was measured: on one `simple_cone` seed the fitter
  returned half-angle 72.85° with its apex 10.8 patch-extents away, explaining the
  six vertices to 6.7e-16 mm, i.e. better than the true cone manages against a
  float32 mesh. Six points genuinely admit more than one cone; Option C escapes
  that only because the ring structure supplies information beyond the six
  positions.
- **Reordering `compute_truncated_cone_axis` to try the exact normal-difference
  method before the apex→centroid guess.** Tried, and much worse — see the DO-NOT
  note on that function. Differencing near-identical normals destroys the signal.
- **Orthogonal-distance geometric fitting** (Lukács et al.; Ahn, *Least Squares
  Orthogonal Distance Fitting of Curves and Surfaces*, Springer LNCS 3151). This
  is the right *refinement* and is essentially what the existing LM already does.
  It does not address initialization, which is where every measured failure is.

### Sources

- Andrews, J.; Séquin, C. H. *Type-Constrained Direct Fitting of Quadric
  Surfaces.* Computer-Aided Design & Applications 10(a), 2013.
- Andrews, J.; Séquin, C. H. *Generalized, basis-independent kinematic surface
  fitting.* Computer-Aided Design 45, 2013, 615–620.
- Lukács, G.; Martin, R.; Marshall, D. *Faithful Least-Squares Fitting of Spheres,
  Cylinders, Cones and Tori for Reliable Segmentation.* ECCV 1998, 671–686.
- Petitjean, S. *A survey of methods for recovering quadrics in triangle meshes.*
  ACM Computing Surveys 34(2), 2002, 211–262.
- Pottmann, H.; Lee, I.-K.; Randrup, T. *Reconstruction of kinematic surfaces from
  scattered data.* 1998.
- Taubin, G. / Pratt, V. — direct algebraic surface fitting, the source of the
  generalized-eigenvalue form used in Options A and B.

---

## Nested fitting: a rank ladder instead of a type ladder (surveyed 2026-08-06)

The question this section answers: **could one fitter, run on a growing patch,
return "the simplest surface that explains what I have so far" — plane, then
cylinder, then cone, then sphere, then torus — with the type an *output* rather
than a precondition?** Asked because that is the missing piece of Phase 1.2. The
earn loop died on a chicken-and-egg — "a cone cannot be fitted from the two
facets a validated loop has when it must first decide, so it never grows" — and
a fitter that never has to decide first does not have that problem.

This is a survey with the option set, the costs, and the two families that are
*not* the answer. **The rank-ladder claim has since been measured** — see "The
rank ladder measured" below, which demotes it from backbone to complement: it is
decisive for sphere and plane and structurally blind to the cylinder. Read the
derivation first; the verdict depends on it.

### The nesting must come from tangency, not from algebra

The obvious nested family is the **general quadric**, and it is the wrong one for
this pipeline, for two reasons that are already measured here rather than
speculative.

It **omits the torus**, the one surface every fillet needs, while including
ellipsoids, paraboloids and hyperboloids that CAD essentially never emits. On
support that barely pins a circular cylinder, the best algebraic quadric is
generically some elliptic cylinder or hyperboloid. Escalating to a quartic to
capture the torus is strictly worse: 34 coefficients, ≥34 points, degree-8
moments, and an even larger zoo of irrelevant surfaces. No reverse-engineering
work in the literature fits tori that way.

And it is **points-only, which is exactly the channel that cannot discriminate**.
The Invariants section measured that a wrong type explains a true support at the
same residual as the right one — a sphere fits 393 of 397 true cylinders *better*
than the cylinder does, at 0.2–0.5 float32 quanta. A points-only algebraic fitter
inherits that blindness whole. Type-projecting the result does not repair it:
a right circular cone requires a repeated eigenvalue in `M`, a condition Andrews
& Séquin note is neither linear nor quadratic in the quadric coefficients, so it
cannot be imposed on the linear fit (this is Option B under "Cone
initialization", with its ≥9 points and its algebraic-distance bias).

What nests correctly is the **space of tangent velocity fields**. Fit the
instantaneous similarity field

```text
v(p) = γ p + c × p + c̄        (7 unknowns: dilation γ, rotation c, translation c̄)
```

and require it tangent to the surface, `v(p)·n(p) = 0`. Each oriented facet
contributes **one linear row**, homogeneous:

```text
γ (p·n) + c·(p × n) + c̄·n = 0
```

This is Pottmann's kinematic surface fitting as generalized by Andrews & Séquin,
both already cited above. What makes it the answer to the question is not the
solve but the **dimension of its null space**, which is a property of the surface
and nothing else:

| surface | tangent similarity fields | null dim | rank | facets to *confirm* |
|---|---|---|---|---|
| plane | 2 in-plane translations, rotation about the normal, one dilation | **4** | 3 | 4 |
| sphere | 3 rotations about the centre | **3** | 4 | 5 |
| cylinder | rotation about + translation along the axis | **2** | 5 | 6 |
| cone | rotation about the axis + **dilation from the apex** | **2** | 5 | 6 |
| torus / any surface of revolution | rotation about the axis | **1** | 6 | 7 |
| extrusion | translation along the generator | **1** | 6 | 7 |
| freeform | none | **0** | 7 | never |

(Confirmation costs one facet more than the rank, for a reason the probe had to
find the hard way — see finding 4 below.)

**The dimension order is the simplicity order**, and it puts the sphere between
plane and cylinder rather than in a separate branch. Under-determination behaves
correctly for free: with `k` rows the rank is at most `k`, so the null space is at
*least* `7 − k`-dimensional, and a patch that has not yet earned its type reads as
a simpler surface than it is. That is precisely "when there's insufficient
support to pin down the result, default to the simplest surface that fits" — as a
**rank statement rather than a residual statement**, which matters because the
residual demonstrably cannot carry it.

Dimension alone does not separate cylinder from cone; the null space's *content*
does, and reading it costs nothing:

| reading | type |
|---|---|
| dim 4 | plane |
| dim 3 | sphere |
| dim 2, every null vector has `γ = 0` | cylinder |
| dim 2, some null vector has `γ ≠ 0` | **cone** — the dilation from the apex |
| dim 1, `c ≠ 0` | torus / surface of revolution |
| dim 1, `c = 0` | extrusion |
| dim 0 | freeform |

Parameters come out of the surviving generators: the rotation generator's Plücker
line **is** the axis, and `A = −c̄/γ` on the cone branch is the apex.

### Why this is cheap here specifically

**Three of the four Phase 1.1 initializers are already sub-blocks of this one
system**, which is the strongest argument that the unification is real rather
than decorative:

- The torus's meridian solve `qᵢ·d + nᵢ·m = 0` (`axis_line_from_normals`) is the
  **rigid sub-block**, `γ = 0`, with `(d, m) = (c, c̄)` read as a Plücker line.
- Cone Option A's `γ(p·n) + a·n = 0` is the **scaling sub-block**, `c = 0`.
- `cylinder()`'s `n₁ × n₂` and `sphere()`'s concurrent circumcenter normal lines
  are 2-row closed forms of the same structure.

The 7×7 Gram matrix *reads which structure is present* instead of assuming one
per solve. It is also **incremental** — one rank-1 update per admitted face, one
7×7 symmetric eigensolve per promotion check — which the per-type LM fitters can
never be, so the check can afford to run at every face admission.

Every measured lesson above transfers unchanged: normalize rows to the patch's own
centroid and extent (the scale-mismatch note under the torus construction — and
the normalization is provably rank-preserving, since translating and scaling the
coordinates maps the 7-parameter space to itself invertibly, so it is conditioning
and never an answer); select the facets to *initialize* from by ruling alignment
and angular spread while *collection* stays validation-gated (the
similarity-vs-spread separation Phase 1.2 paid for); and treat the result as a
**seed for the vertex-based least squares**, never as the answer, because every
row is built from a facet normal and inherits the `O(Δθ²/8)` chord-plane lean.

**Why the facet centroid is the right evaluation point**, and it is not obvious:
the row wants a point on the surface paired with the true normal, and the facet
gives neither exactly. For a ruling-aligned quad on a cylinder the two errors
cancel — the chord of a circular arc is parallel to the tangent at the arc's
*midpoint*, and the centroid sits at that midpoint angle, so the rotation
generator annihilates the row exactly. Evaluating at a vertex instead is `O(Δθ)`
wrong, an order worse. On a cone the facet plane contains the apex exactly, so the
`p·n` entry is exactly zero however the facet is shaped. On a **sphere** the exact
point is the facet's **circumcenter** rather than its centroid, since that is
where the normal line passes through the centre — the same distinction Phase 1.1
recorded for `sphere()`. Both are measured under family K below.

### What it does not fix

Stated plainly, because the failure modes are known and the rank ladder does not
touch them:

- **Genuine ambiguity stays ambiguous.** A single strip band — vertices on two
  circles, every facet centroid at one axial position — is sphere-compatible in
  the normal channel as well as the vertex channel (`simple_cylinder`'s 32-quad
  band, residual 2.07e-7 mm). The rank signal is degenerate there for the same
  reason the residual is. *(Measured, and far worse than "stays ambiguous": the
  band's rank collapses all the way to a plane's, and it is the common case
  rather than a corner one — see finding 1 below.)* The load stays on the structural gates
  (`MIN_AXIAL_SWEEP_DEG`, solid-angle coverage) and on branch-don't-choose.
- **Rank is not to be read with a threshold.** The eigenvalue spectrum orders
  candidate interpretations; **validation decides**. This is the same rule
  `MIN_NULLSPACE_SEPARATION` is documented under, arrived at from the other end —
  the measured sphere-vs-torus conditioning gap is 2.06e-6 against 9.2e-4, a
  factor of 450 with pure data-dependence in between, so no cut exists.
- **Support cost is real and is the headline risk.** A cylinder needs 6
  independent rows against the 2 facets `initialize::cylinder` needs, a torus 7
  against a vertex star's 4–12. The rank ladder is a *classifier* that needs more
  support than the closed forms need, and whether the extra support is available
  before a patch crosses a junction is the whole question family K asks.
  *(Measured: it is worse than this, because the cylinder also needs its support
  spread over two axial bands, which a ruled tessellation does not supply.)*
- **Plane starvation of the freeform rung.** An earned plane patch is by
  construction "the faces that validate a plane", so a single-facet plane always
  succeeds — measured to collapse the NURBS models to 94/51/83 per-facet planes.
  Whatever consumes the ladder must keep the leftover-coherence policy.

### The option set, ranked (before measurement; see the verdict below)

- **A — the kinematic rank ladder above.** Proposed as the backbone. Adds the
  normal channel the points-only methods lack, subsumes three existing
  initializers, incremental, and its degeneracies are the ones already
  understood. **Measured, and demoted — see "The rank ladder measured" below.**
- **B — general quadric with eigen-projection.** Demoted to triage
  (analytic-vs-freeform) and initializer of last resort, for the reasons above.
- **C — per-type fitters with degeneracy-safe parameterizations.** The
  conservative option, and *complementary rather than alternative*: adopt
  Lukács–Martin–Marshall's parameterization discipline, whose stated property is
  that "as the principal curvatures decrease, the results naturally become closer
  to the surfaces of simpler type, unlike other methods which may diverge". That
  is the same move as this document's cone `taper = tan α` and the bent tube's
  `k = 1/R`. Keep it in the refiners whatever backbone wins; it leaves the four
  initializers un-unified and promotion hand-enumerated if used alone.
- **D — Dupin cyclides.** The one classical family containing *exactly*
  {plane, sphere, cylinder, cone, torus} and closed under the relevant
  degenerations, with a genuine unified theory in the space of oriented spheres.
  Declined: the literature is construction- and blend-oriented, mesh fitting is
  thin, the Minkowski–Lorentz machinery is heavy, and the extra coverage is a
  liability in a pipeline whose whole discipline is "simplest surface that
  validates". Worth remembering only if variable-radius blends become a target.

### The rank ladder measured (family K)

`src/stage2/kinematic.rs` builds the system; `--invariant-report` family K scores
it against `--compare` ground truth. Per reference face, over a BFS patch grown
from several starts and sampled at a ladder of support sizes: the 7 relative
eigenvalues, the null dimension read from the widest ratio gap *that row
starvation does not already explain*, the type that dimension plus its
null-space content names, whether that matches the truth, the straddling
fraction and axial band count of the support, and — the analogue of family F —
the residual
of the **reference's own** tangent generators against the measured rows, at both
the centroid and the circumcenter, which asks whether the signal is there to be
read before asking whether the read succeeds.

```bash
./tools/invariants.sh --out tmp/kin.tsv --jobs 8 --kinematic && python3 tools/kin_analyze.py tmp/kin.tsv
```

**Measured 2026-08-06 over 138 models, 1360 reference faces, 62 774 readings.**
Growth confined to each reference face's own facets, so this is the signal with
*no* junction contamination — the best case, deliberately.

Headline, per reference face: does the ladder read the true type at *any*
support size on the ladder 3…48?

| true type | faces | best evaluation point | reads correctly | first size that works (p50) |
|---|---|---|---|---|
| plane | 374 | either | **94.9%** | 4 |
| sphere | 133 | **circumcenter** | **100.0%** | 5 |
| cone | 83 | **centroid** | 66.3% | 6 |
| torus | 131 | either | 71.0% | 12 |
| **cylinder** | 639 | — | **4.5%** | — |

**The verdict: the ladder is not a backbone, because it is structurally blind to
the cylinder — and the cylinder is 47% of the corpus's analytic faces.** The
other four types work, sphere and plane decisively. The remaining findings say
why, and two of them are worth more than the headline.

**1. The cylinder failure is structural, and it is the tessellation's doing.**
In a band of ruling-aligned quads every facet centroid sits at the same axial
position `z_c` and the same radial distance `r_c`, so `p·n = r_c` is one constant
and `p × n = z_c (d × n)` is a fixed multiple of the normal: **every row is an
affine function of `n` alone, the row space is spanned by {constant, nₓ, n_y},
and the rank is 3 — a plane's rank — at any circumferential extent.** Measured
rather than inferred, by counting distinct axial centroid positions alongside the
rank:

| distinct axial bands in the support | readings | mean rank | rank = 3 | read correct |
|---|---|---|---|---|
| **1** | 2646 | 3.35 | **64.9%** | **0.0%** |
| 2 | 349 | 6.00 | 0.0% | 5.7% |
| 3+ | 2567 | 6.18 | 3.4% | 0.7% |

and a cylinder is *ruled*, so a CAD tessellator emits exactly one row of quads
spanning its whole height — the same fact this document already records for the
cone ("one row of quads spans its entire height"). 56–61% of the ladder's wrong
cylinder readings are literally `plane`. This is the axial analogue of "a short
arc cannot determine a free surface of revolution": **a single band carries no
axial information at all**, and no conditioning fixes it.

The failure is also *diagnosable with no ground truth* — count distinct axial
centroid positions — and the case it fails on is exactly the case
`initialize::cylinder`'s `n₁ × n₂` handles best, because one band of quads is
100% ruling-aligned. Measured on the same faces in the same run: the closed form
lands the axis under **0.1° on 89.7%** of 261 cylinder reference faces, p50
0.0000°. **They are complementary, and the complementarity is exact rather than
lucky.**

**2. The evaluation point is type-dependent, decisively, and the derivation was
right.** The reference's own generators against the measured rows (relative
residual, p50):

| true type | centroid | circumcenter | reads right: centroid → circumcenter |
|---|---|---|---|
| plane | 0.00 | 0.00 | 90.6% → 94.7% |
| **sphere** | 3.26e-3 | **3.88e-7** | 83.5% → **100.0%** |
| cylinder | 1.40e-6 | 1.43e-6 | 2.5% → 4.4% |
| **cone** | 4.91e-4 | 1.29e-4 | **62.7% → 32.5%** |
| torus | 4.11e-7 | 3.53e-7 | 64.9% → 67.2% |

Sphere wants the circumcenter by four orders of magnitude, and the cone wants the
**centroid** — note that the cone's generator residual is *better* at the
circumcenter while its *read* is twice as good at the centroid, so residual and
readability are not the same quantity here either. A ladder that fixed one
evaluation point would be measuring its own choice.

**3. The gap measures how cleanly a rank is determined, not whether it is the
right rank.** Split by whether the read was correct: **for every type but the
plane the wrong reads' p95 gap exceeds the right reads' p50**, so no cut
separates them — and on three of the seven type/evaluation-point combinations the
wrong reads' *median* gap is the larger one outright (sphere at the circumcenter,
right p50 2.28e9 against wrong 4.73e10; cone at the centroid, 8.35e3 against
9.05e6; cone at the circumcenter, 5.80e7 against 4.21e8). A one-band cylinder is
confidently, enormously wrong — rank 3 with a gap of ~1e12. So the standing rule
against thresholding a spectrum arrives from a third direction, and in its
sharpest form yet: **a large gap is evidence that a rank is well determined and
no evidence at all that it is the right rank.**

**4. Two defects the probe found in its own reading rule**, both recorded because
both are the obvious implementation:

- *Reading the widest gap over the whole spectrum measures `7 − rows`.* With `k`
  rows the Gram matrix has `7 − k` zero eigenvalues for reasons having nothing to
  do with the surface. The first version walked 4, 3, 2, 1, 0 as support went 3,
  4, 5, 6, 8 and hit the true type exactly once per model, by coincidence, where
  `7 − k` happened to equal it. The fix is to compare the read against that
  trivial floor, which costs one row more than the type's rank: 4 to confirm a
  plane, 5 a sphere, 6 a cylinder or cone, 7 a torus — and freeform never
  confirmable, only "nothing simpler was confirmed", which is its correct status.
- *A content threshold sharp enough to name a type is wrong most of the time.*
  `BLOCK_PRESENT` was documented as safe because "a cylinder's `γ` is identically
  zero and a cone's is of order one, so any cut well below 1 gives the same
  answer". On real meshes the rank reads a true torus's dimension correctly
  **79–81%** of the time while that floor labels **70%** of those readings
  `Spiral`, because a facet normal's `O(Δθ²/8)` error lands ~1e-3 in the `γ`
  block. `Spectrum::candidates` now returns every type a dimension admits,
  ordered, and validation picks. Separating the two numbers is what makes the
  torus row above readable at all — note this table is **per reading** at its
  best evaluation point, where the headline table is per *face* over all support
  sizes, so the two answer different questions and the torus's 81.5% and 71.0%
  are not in conflict:

| true type | dimension read correctly | top-1 type correct |
|---|---|---|
| plane | 95.9% | 95.9% |
| sphere (circumcenter) | 97.4% | 97.4% |
| cone (centroid) | 28.6% | 28.6% |
| **torus** | **81.5%** | 11.2% |
| cylinder | 5.0% | 0.6% |

**5. Failure direction is mostly right, but not entirely.** The design requires
that an unearned read degrade to something *simpler*, never to a different type
at equal complexity. Of the wrong cylinder readings, 57–61% are simpler and
26–28% are under-determined — but **12–16% name a wrong type outright** (mostly
sphere), and for the torus that figure is 47–49%. So "under-determination
degrades to simplicity" holds as a tendency, not as a guarantee, and anything
built on the ladder still needs validation to catch the rest.

### What this changes

- **Option A is demoted from backbone to complement.** Use the rank ladder where
  it is decisive — **sphere** (100.0%, circumcenter, from 5 facets) and **plane**
  (94.9%, from 4) — and keep `initialize::cylinder` for the cylinder, where it is
  90% accurate on exactly the configuration the ladder cannot read. The torus's
  81.5% dimension read is a usable axis signal and is close to what
  `axis_line_from_normals` already extracts, which is unsurprising: that solve is
  this system's rigid sub-block.
- **The cone stays the weak type**, at 66.3% here against the 49.4% of the
  Invariants section's seed rules — an improvement, not a solution, and one that
  needs the centroid evaluation to get even that.
- **The support cost is the real barrier to using this inside a growth loop.**
  Confirming a cylinder or cone needs 6 rows and a torus 7, against the 2 facets
  `initialize::cylinder` needs and the 4–12 of a torus vertex star. A patch that
  must reach 7 validated facets before it can name its own type is a patch that
  may already have crossed a junction — which is the Phase 1.2 chicken-and-egg
  arriving in a new costume rather than being dissolved.

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

**Centring matters, and its absence has a hard failure point (measured
2026-07-31).** For a *cone* the uncentred covariance is not safe. A cone's
normals satisfy $n = \cos\alpha\,\hat r - \sin\alpha\,d$, so over a full
revolution $M$ has eigenvalues $\cos^2\!\alpha/2$ (twice, perpendicular to the
axis) and $\sin^2\!\alpha$ (once, along it). The axis is the *smallest*
eigenvector only while $\sin^2\!\alpha < \cos^2\!\alpha/2$, i.e.
**$\alpha < 35.26°$** — above that it becomes the *largest* and the estimate
silently flips. This is not exotic: a 90° countersink is $\alpha = 45°$, and the
corpus's own `cone_15x20` is 36.87°.

**Centring removes the regime entirely.** Subtracting the mean normal cancels
the $\sin\alpha\,d$ term and leaves $\cos\alpha(\hat r_i - \bar{\hat r})$,
which lies in the plane perpendicular to the axis **for every half-angle**, so
the smallest eigenvector of the centred covariance is the axis exactly. Over a
*partial* arc the perpendicular spread is also anisotropic, which degrades the
uncentred estimate well before 35.26°; centring is unaffected. Use the centred
form.

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
- **Cone**: eigenvalues $(\sin^2\theta, -\cos^2\theta, -\cos^2\theta)$ — one of opposite sign to a repeated pair, so $\det(\mathbf{A}) \neq 0$. What is singular for a cone is the *extended* $4 \times 4$ matrix, because the apex is a singular point of the surface. (An earlier version of this line said "three nonzero eigenvalues with $\det(\mathbf{A}) \approx 0$", which is the **cylinder** condition — a cylinder is the case with one zero eigenvalue.)
- **Plane**: $\mathbf{A} = 0$ (all eigenvalues zero).

Quadric fitting is fast (SVD of an $N \times 10$ matrix, $O(N)$) and provides a good initialization for specialized fitters. It could be used as an alternative initialization for cone fitting — see **Option B** under "Cone initialization: the options not yet built" above for what that costs (≥9 points, algebraic-distance bias) and for the two alternatives.
