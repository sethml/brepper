# Surface Fitting Algorithms

This document describes the surface fitting algorithms used in brepper's Stage 2 pipeline. Each section covers the hypothesis data structure, seeding strategy, BFS region growing, fitting math, and validation for one surface type.

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for the overall pipeline architecture and implementation status.

---

## General Approach

"Region growing" from seed faces while preserving fit invariants, refitting as necessary as faces are added to the region. Given that vertices lie within epsilon of precisely on the surfaces involved, any seed on a surface should grow to the entire region.

---

## 2.1 Deduce planar hypotheses
Fit planar hypotheses to all connected sets of coplanar faces. Every face is assigned to exactly one planar hypothesis (including single-face hypotheses for faces on curved surfaces). A planar hypothesis consists of:
- `normal: [f64; 3]` — Unit normal vector pointing outward from the shell/solid.
- `distance: f64` — Signed distance from origin to plane along the normal (plane equation: normal · p = distance).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this plane (unordered set, not right-hand-rule ordered — boundary representation is deferred to stage 3).
- `error_max, error_min: f64` — Maximum (positive) and minimum (most negative) signed distance from any vertex to the plane.
- `error_abs_sum: f64` — Sum of absolute vertex-to-plane distances.

The algorithm uses BFS region growing (not DFS as originally proposed — BFS explores more uniformly and avoids the "narrow corridor drift" problem where DFS could accumulate error along a thin strip before reaching the main planar region):
- Initialize all `planar_hypothesis` to `UNDEDUCED_PLANAR_HYPOTHESIS` (-2).
- For every face index fi:
    - If `face[fi].planar_hypothesis != -2`: continue (already assigned).
    - Create a new planar hypothesis hi. Seed it with face fi's normal and compute the initial plane distance as the average dot product of fi's vertices with the normal.
    - Assign fi to hypothesis hi. Add fi to a BFS queue.
    - BFS loop — pop faces from the queue and examine each neighbor ni:
        - If ni is already assigned to a hypothesis, skip.
        - **Vertex distance check**: compute signed distance from each vertex of ni to the current plane. If all are within `--vertex-tolerance` (default 1e-5 mm), accept directly. If any vertex exceeds `2 * vertex_tolerance`, skip immediately — the face is clearly non-coplanar and re-fitting cannot help (since the plane can shift by at most ~tolerance while keeping existing vertices in tolerance). (No separate angular alignment check is needed — with tight vertex tolerance, any triangle with all vertices in the tolerance band necessarily has a closely aligned normal.)
        - **Re-fit attempt**: if any vertex exceeds tolerance but none exceeds `2 * tolerance`, try re-fitting the plane using area-weighted normal averaging over all current faces plus ni, with the distance as the average projection of all current+new vertices. Then check that **all** vertices (existing and new) are within tolerance of the re-fitted plane. If so, accept the re-fit; otherwise skip ni. (This addresses the comment's concern about previously-accepted faces no longer fitting — the re-fit validation checks all vertices, not just the new face's.)
        - If accepted: assign ni to hypothesis hi, add its vertices to the vertex set, push ni onto the BFS queue.
    - After BFS completes: final re-fit using all collected faces and vertices (area-weighted normal averaging + vertex centroid). Compute error metrics (error_max, error_min, error_abs_sum).
    - Single-face hypotheses are kept (not deleted). They represent faces on curved surfaces that will be assigned to cylindrical/spherical hypotheses in later stages.

Plane fitting uses **area-weighted normal averaging**: each face's unit normal is weighted by the face's triangle area, then the weighted sum is normalized. The plane distance is the mean of `normal · vertex` over all vertices in the hypothesis. This gives larger faces more influence on the plane orientation, which is more robust than unweighted averaging.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted plane (perpendicular projection).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - (Changed from original plan: the plan described projecting face *edges* onto the hypothesis and checking the enclosed surface region. The centroid-projection approach is simpler and sufficient for validation — it confirms the fitted plane aligns with a real STEP surface at representative interior points. Edge projection would be needed for boundary accuracy validation, which is a stage 3 concern.)

*Comment on 2.1: ~~The DFS vs BFS question is worth considering: DFS can get "trapped" in a narrow corridor and accumulate drift before exploring the main region. BFS explores more uniformly. However, the re-fitting step partially mitigates this.~~ Resolved: BFS was chosen for the implementation. ~~A bigger issue: once you re-fit the plane, previously accepted faces might no longer fit the new plane! Consider a final validation pass that removes faces whose vertices exceed the tolerance after the final fit.~~ Resolved: the re-fit step checks all vertices (existing + new) before accepting, and a final re-fit with error metrics is computed at the end. ~~Also: the "vertices in right-hand-rule order" for the hypothesis is unclear—planar regions aren't simply connected in general (they can have holes), so you'll need a more complex boundary representation.~~ Resolved: vertices are stored as an unordered set; boundary representation is deferred to stage 3.*

---

## 2.2 Deduce cylindrical hypotheses
Fit cylindrical hypotheses to connected sets of faces that lie on a common cylinder. A cylindrical hypothesis consists of:
- `axis_origin: [f64; 3]` — A point on the cylinder axis.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis.
- `radius: f64` — Radius of the cylinder (always positive).
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a pipe) or toward it (concave, like the inside of a hole).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this cylinder.
- `error_max: f64` — Maximum absolute distance from any vertex to the cylinder surface (i.e. max of `| ||v - axis_closest|| - radius |`).
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stage 2.1 but seeded from triples of nearby faces (since a single triangle is always planar and cannot determine curvature, and two faces provide minimal angular constraint on the cylinder axis). A key challenge is that **any two non-coplanar faces define some cylinder** — the cross product of their normals gives an axis, and a circle fit always produces a radius. The previous algorithm took the first valid seed and committed to it, which meant a bogus 3-face cylinder from an unfortunate seed pair could permanently consume faces that should belong to a correct 20-face cylinder seeded from a different pair. The fix is multi-seed evaluation with 3-face seeds that span more of the cylinder's circumference:

**Pairwise qualification:**
Two faces fi and ni are pairwise qualified to participate in a cylindrical hypothesis if:
- fi and ni are neighbors.
- fi and ni are unassigned (`cylindrical_hypothesis == -2`).
- fi and ni have a sufficiently different normal (cross product magnitude exceeds `MIN_CROSS_THRESHOLD`, e.g. 0.01).
- The dihedral angle between fi and ni does not exceed `--angular-tolerance` (default 17.5°).

**Cylinder fit qualification:**
A face or set of faces is cylinder fit qualified if:
- All vertices are within `--vertex-tolerance` of the cylinder.
- All face centroids are within `--surface-tolerance` of the cylinder.
- **Convexity check**: compute the vector from the axis to the face's centroid. Check that the face's normal agrees with the hypothesis convexity (dot product with radial vector has the expected sign). Reject if inconsistent — this prevents merging the inner and outer surfaces of a thin-walled cylinder.

**Multi-seed evaluation:**
- Initialize all `cylindrical_hypothesis` to `UNDEDUCED` (-2).
- For each face fi where `cylindrical_hypothesis == -2`:
    - For each neighbor n1 of fi where n1 is pairwise qualified with fi:
        - For each neighbor n2 of n1 where n2 is not a neighbor of fi and n2 is pairwise qualified with n1:
            - Produce a cylinder fit to the set of vertices belonging to fi, n1, and n2.
            - If faces fi, n1, and n2 are cylinder fit qualified, accept this set of faces as a valid seed and perform BFS expansion on them using temporary data structures (a local face list and vertex set — do not mutate `mesh.faces[].cylindrical_hypothesis` during trials). The initial BFS queue should contain fi, n1, and n2, so we explore all of their neighbors. Each trial produces a candidate hypothesis with a face count and fitted parameters. After each trial completes with acceptance, compare against the **current best candidate** for this starting face; keep whichever has more area covered. This avoids storing all trial results.
    - If no valid (fi, n1, n2) triple was found, skip fi (leave as `UNDEDUCED` — it may be absorbed by a later face's BFS, or will remain unassigned at the end).
- Accept the **best candidate** as a cylindrical hypothesis: create a cylindrical hypothesis, and for every face in the hypothesis, assign the face's cylindrical hypothesis field.
- Apply **angular coverage validation** (see below) to each trial. If it fails, discard the trial.
- Accept the **best candidate** as a cylindrical hypothesis: create a cylindrical hypothesis, and for every face in the hypothesis, assign the face's cylindrical hypothesis field.
- After all faces have been processed, set any remaining `UNDEDUCED` faces to `NO_HYPOTHESIS`.

Future work (not currently implemented - in the future we'll evaluate whether these optimizations help):
    - **Early termination** of redundant trials: during a trial BFS, if a best candidate already exists, check whether the trial is rediscovering the same cylinder. After the trial adds a face, compare its fitted cylinder parameters against the best candidate's: axis directions nearly parallel (`|a_trial · a_best| > 1 - 1e-6`), radii within `--vertex-tolerance`, and axis-to-axis distance within `--vertex-tolerance`. If parameters match AND every face accepted so far is already in the best candidate's face set, abandon the trial — it can only produce a subset of the existing best. This is the dominant optimization: on a 20-face fillet cylinder with ~3 valid seeds per face, most of the ~60 trials terminate after 2–3 faces instead of exploring all 20.

**BFS expansion** (run once per trial seed):
- Repeatedly pop a face fi from the queue and examine each neighbor ni:
    - If ni is already assigned a cylindrical hypothesis (in the mesh) or claimed by this trial, skip.
    - **Angular tolerance check**: for each of ni's mesh neighbors that is already assigned in this trial, compute the dihedral angle between ni and that neighbor. If any exceeds `--angular-tolerance`, skip ni. Checking all assigned neighbors (not just the BFS parent) provides defense-in-depth against creased NURBS surfaces where BFS might approach a face from a low-angle direction while a high-angle assigned neighbor exists on a different axis. (Note: the full pairwise qualification — including the cross product minimum and unassigned requirement — is NOT applied during BFS expansion. The cross product minimum is a seeding criterion only: on a finely tessellated cylinder, adjacent faces legitimately have near-parallel normals with tiny cross products. The unassigned criterion is N/A since trial-claimed faces are logically assigned within the trial.)
    - If ni fails cylinder fit qualification:
        - If any vertex of ni is further than `2 * vertex_tolerance` from the cylinder or ni's centroid is further than `2 * surface_tolerance` from the cylinder, skip ni. Otherwise, do a re-fit attempt:
        - **Re-fit attempt**: re-fit the cylinder from all current faces plus ni (see Cylinder Fitting below). Check that **all** faces (existing and new) are cylinder fit qualified with the re-fit. If so, accept the re-fit; otherwise reject the re-fit cylinder parameters and skip ni.
    - If accepted: add ni to the trial's face list and vertex set, push ni onto the BFS queue.
- After BFS completes: do a final re-fit from all accumulated faces and vertices. Compute error metrics.

**Angular coverage validation** (applied to each trial):

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

**Computational cost**: For each starting face with K valid seed partners (neighbors passing pairwise qualification), each partner n1 may have L additional neighbors qualifying as n2, producing up to K×L trial BFS runs per face. Typical K is 1–3 (most faces have 3–4 neighbors, of which 1–2 have sufficiently different normals) and L is 1–2 (n1's neighbors minus fi and minus n1's non-qualifying neighbors). So typically 1–6 trials per face, each of which is a BFS that covers the connected cylinder region. Without early termination, total work is proportional to `sum over faces of (K_i × L_i × region_size_i)`. Since most faces on the same cylinder will discover the same region, the work per cylinder is roughly `N_faces × K × L × N_faces`, dominated by the first few faces that discover the region before their neighbors get committed. The temporary data structures (HashSet for claimed faces, Vec for face list) are lightweight.

**Cylinder fitting** (used for seeding and re-fitting):

The fitting proceeds in three steps:

1. **Axis direction estimation from face normals.** On a cylinder, every face normal is perpendicular to the axis: `n · a = 0`. So the axis direction `a` minimizes $\sum w_i (n_i \cdot a)^2$ subject to $\|a\| = 1$, where $w_i$ is the area of face $i$. This is the eigenvector corresponding to the *smallest* eigenvalue of the 3×3 weighted covariance matrix $M = \sum w_i \, n_i \, n_i^T$. For the three-face seed, this is the full 3×3 eigenvector computation (the cross-product shortcut `a = normalize(n1 × n2)` that works for two-face seeds is not used, since we want the axis estimate that best fits all three faces' normals).

2. **Axis position and radius via 2D circle fitting.** Given the axis direction `a`, choose two orthogonal unit vectors `u`, `w` perpendicular to `a`. Project each vertex onto the 2D plane: `(x_i, y_i) = (v_i · u, v_i · w)`. **Center the 2D coordinates** by subtracting the mean: `x_i' = x_i - mean_x`, `y_i' = y_i - mean_y`. Fit a circle to the centered 2D points using the algebraic least-squares method: fit the linear model $x^2 + y^2 + Dx + Ey + F = 0$ via least squares, then `center' = (-D/2, -E/2)` and `radius = sqrt(center_x'² + center_y'² - F)`. The axis origin in 3D is `(center_x' + mean_x) * u + (center_y' + mean_y) * w` (the component along `a` is arbitrary). Centering is essential for numerical stability: without it, the normal equations matrix becomes ill-conditioned when absolute vertex coordinates are large relative to the arc span (e.g., vertices at (15mm, 12.5mm) spanning 2.5° of a r=2mm cylinder produce sums-of-squares ~1350 that dwarf the r=2 signal, causing wildly wrong radii).
3. **Levenberg-Marquardt refinement.** The initial axis + circle fit provides a good starting point but can have ~0.4° axis error when faces span a narrow arc (<15°). LM refines all 6 parameters simultaneously: `[alpha, beta, qx, qy, qz, radius]` where alpha/beta are tilt angles from the initial axis in two perpendicular directions. The residual for each vertex is `||cross(p - q, a)|| - r` (signed distance to cylinder surface). Numerical Jacobian via central differences (eps=1e-8). Uses `levenberg-marquardt` crate v0.14 with `nalgebra` v0.33 types. Falls back to the step-2 estimate if LM produces a negative radius.

The normal covariance matrix `M` can be maintained incrementally during BFS (add $w_i \, n_i \, n_i^T$ when accepting a face), but the eigenvector solve, circle fit, and LM refinement must be recomputed from scratch when re-fitting, since a change in axis direction invalidates the 2D projection. Steps 1-2 are O(n) in vertices with simple arithmetic; step 3 adds ~5-10 LM iterations (typically <1ms per hypothesis).

**Subtleties and edge cases:**
- **Near-parallel seed normals**: Pairwise qualification requires cross product magnitude > `MIN_CROSS_THRESHOLD` (e.g. 0.01, ~0.6° between normals) for each consecutive pair in the seed chain (fi↔n1 and n1↔n2). This ensures each pair contributes meaningful angular information. Note that fi and n2 need NOT have different normals — on a cylinder, faces at the same circumferential position but different axial positions have identical normals. The angular diversity comes from the intermediate face n1 having a different circumferential position.
- **Cones vs cylinders**: A cone has varying radius along the axis. If we fit a cylinder to conical faces, vertex distances near the ends will exceed tolerance, so cone faces will be naturally rejected. Cone fitting is deferred.
- **Large-radius cylinders**: A cylinder with very large radius looks locally planar. Such faces will belong to multi-face planar hypotheses from stage 2.1 (since vertices on a nearly-flat cylinder are nearly coplanar within tolerance). Multi-face planar faces CAN seed cylinders — the angular coverage validation, minimum face count, and other checks prevent bogus large-radius cylinder fits. Stage 2.6 resolves any ambiguity between cylindrical and planar assignments.
- **Cylinder wrapping >180°**: This works naturally with BFS. The convexity check ensures we don't accidentally merge inner and outer surfaces.
- **Multiple cylinders sharing an axis**: E.g. stepped bore holes with different radii. These will form separate hypotheses because vertices at different radii will exceed the vertex tolerance.
- **Seam between cylinder and plane**: Faces at the junction have some vertices on the plane and some on the cylinder. Since the plane vertices aren't on the cylinder surface, the vertex distance check rejects them. This correctly excludes boundary faces.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the cylinder surface (nearest point: project onto axis, then move radially to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta (the flat triangle's centroid lies inside/outside the curved surface by approximately $h^2 / 8R$ where $h$ is the chord length). For typical tessellations this is small relative to `--surface-tolerance`, but it could become significant for very coarse meshes on tight-radius cylinders.

---

## 2.3 Deduce spherical hypotheses

Fit spherical hypotheses to connected sets of faces that lie on a common sphere. Faces with cylindrical hypotheses are not excluded — a face can legitimately belong to both a cylinder and a sphere (e.g., equator faces). Stage 2.6 resolves overlapping hypotheses later. A spherical hypothesis consists of:
- `center: [f64; 3]` — Center of the sphere.
- `radius: f64` — Radius of the sphere (always positive).
- `convex: bool` — Whether face normals point away from the center (convex, like the outside of a ball) or toward it (concave, like a bowl).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this sphere.
- `error_max: f64` — Maximum absolute distance from any vertex to the sphere surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stages 2.1 and 2.2 but seeded from the group of faces surrounding a mesh vertex. A vertex neighborhood on a sphere naturally provides normals spanning 3D (unlike cylinder seeding which only needs normals spanning a 2D plane), making it a more natural seed geometry for the 4-DOF sphere fit. The algebraic sphere fit from 4+ non-coplanar vertices determines a unique sphere, and validation catches bad fits:

**Seeding:**
- For each mesh vertex V, collect the faces incident on V ("surrounding faces"). Exclude any face that has `spherical_hypothesis != UNDEDUCED` or that belongs to a multi-face planar hypothesis (same rationale as cylindrical seeding: their vertices lie on flat surfaces, not spheres). Multi-face planar faces remain UNDEDUCED so BFS can absorb them later.
    - If fewer than 3 non-excluded faces remain around V, skip this vertex (not enough angular diversity for a reliable sphere fit).
    - Compute a spherical fit to all vertices of the remaining surrounding faces using least-squares sphere fitting (see Sphere Fitting below).
    - Skip this vertex if:
        - The fitted radius exceeds `max_sphere_radius` (see below).
        - Any vertex of the seed faces is not within `--vertex-tolerance` of the estimated sphere surface.
        - The dihedral angle between any pair of adjacent seed faces (i.e., seed faces sharing a mesh edge) exceeds `--angular-tolerance` (default 17.5°).
        - The distance from any seed face centroid to the sphere exceeds `--surface-tolerance`.
    - Determine convexity: for each seed face, compute the dot product of the face normal with the vector from sphere center to face centroid. Positive → convex, negative → concave. If convexity is inconsistent across seed faces, skip this vertex.
    - Create a new spherical hypothesis, assign all seed faces, and perform BFS expansion.

*Design note: The original plan called for pair-based seeding (two adjacent faces with non-parallel normals, as in cylindrical), then 3-face seeds with a normals-span-3D eigenvalue check. Pair-based seeding failed on finely tessellated spheres where adjacent faces have nearly identical normals (spanning only ~6°). Vertex-neighborhood seeding is superior for spheres: it naturally collects faces spanning multiple directions around a vertex, providing better constraints for the 4-DOF sphere fit without needing an explicit normals-span check. The sphere fit itself validates geometry — if the vertices don't lie on a sphere, the fit produces a large radius or high error, caught by max_sphere_radius and vertex distance checks.*

**BFS expansion:**
- Pop faces from the BFS queue and examine each neighbor ni:
    - If ni already has a spherical hypothesis, skip.
    - **Vertex distance check**: for each vertex of ni, compute `| ||v - center|| - radius |`. If all are within `--vertex-tolerance`, proceed to remaining checks (no re-fit needed). If any exceeds `2 * vertex_tolerance`, skip (re-fitting cannot help, same reasoning as planar and cylindrical). If between 1x and 2x, flag for re-fit.
    - **Centroid validation**: compute the distance from ni's centroid to the fitted sphere. If beyond `2 * surface_tolerance`, skip. If beyond `surface_tolerance`, flag for re-fit. (Same structure as stage 2.2 cylindrical BFS.)
    - **Convexity check**: verify ni's normal agrees with the hypothesis convexity (dot product of normal with radial vector has the expected sign). If not, skip.
    - **Angular tolerance check**: for each of ni's mesh neighbors that is already assigned to this hypothesis, compute the dihedral angle between ni and that neighbor. If any exceeds `--angular-tolerance`, skip.
    - **Re-fit attempt**: if vertex distance or centroid validation flagged a re-fit, re-fit sphere from all current faces plus ni. If all vertices are within `--vertex-tolerance` and all face centroids within `--surface-tolerance` and radius within `max_sphere_radius`, accept the face with the updated sphere parameters. If not, revert to the previous sphere fit and skip the face.
    - If accepted: assign ni to the spherical hypothesis, and add its neighbors to the BFS queue.
- After BFS completes: final re-fit, compute error metrics.
- **Minimum face count**: require at least 4 faces (a sphere needs at least 4 non-coplanar points for a unique fit).

**Solid-angle coverage validation** (applied to the completed hypothesis, analogous to angular coverage for cylinders):

The problem: sphere BFS can grow along cylinder fillet surfaces, because locally adjacent cylinder faces fit on a sphere within vertex_tolerance. This produces oversized sphere hypotheses spanning a thin strip around a cylinder fillet (observed on `onshape_rounded_cube`). The analog of the cylindrical angular coverage test is a **solid-angle coverage test**: faces on a genuine spherical patch produce centroid-to-center directions that span 3D, while faces on a cylinder fillet produce directions that span only a 1D strip.

Algorithm:
1. For each face in the hypothesis, compute the unit vector from sphere center to face centroid: $d_i = \text{normalize}(\text{centroid}_i - \text{center})$.
2. Compute the area-weighted 3×3 covariance matrix of these unit vectors: $C = \sum w_i \, d_i \, d_i^T$ where $w_i$ is the face's mesh area.
3. Compute eigenvalues of $C$. Let $\lambda_1 \geq \lambda_2 \geq \lambda_3$.
4. Require $\lambda_3 / \lambda_1$ to exceed a minimum threshold (e.g., 0.01). If it fails, the hypothesis directions are nearly coplanar (a strip on the sphere, characteristic of cylinder fillet growth) or nearly collinear (a narrow band), and the hypothesis should be discarded.

*Note: Fillet faces lie along a 1D arc on the sphere's direction space, producing $\lambda_3 \approx 0$. Genuine spherical patches (hemispheres, domes, spherical pockets) have centroid directions spanning 3D, so all eigenvalues are substantial. For partial spherical caps, the directions span at least 2D; very small caps may fail but such small patches are unlikely to be meaningful spherical features.*

*Note: Like cylindrical angular coverage, this should NOT be applied during BFS growth — early in BFS, faces are clustered near the seed and naturally appear low-dimensional. It is a global structural validation for completed hypotheses.*

**Maximum radius**: `max_sphere_radius`. With solid-angle coverage validation and surface-tolerance validation during BFS, the radius limit no longer needs to be tight — those checks prevent the pathological growth patterns that the original 10× bounding box diagonal limit was designed to catch. A large value like `bounding_box_diagonal * 1000` suffices as a numerical guard rail, preventing sphere fits with absurd radii (where floating-point precision of vertex coordinates could make flat faces appear spherical) while not rejecting any plausible real-world spherical feature. Compute the bounding box once at stage 2 entry.

**Sphere fitting** (used for seeding and re-fitting):
- Given vertices on a sphere: $|v - c|^2 = r^2$, expand to $v \cdot v - 2 v \cdot c + |c|^2 = r^2$.
- Rearrange: $v \cdot v = 2 v \cdot c - |c|^2 + r^2 = 2 v \cdot c + k$ where $k = r^2 - |c|^2$.
- This is linear in unknowns $(c_x, c_y, c_z, k)$. Solve via least squares: $A x = b$ where $A_{i} = [2v_{ix}, 2v_{iy}, 2v_{iz}, 1]$ and $b_i = v_i \cdot v_i$.
- After solving: center $= (c_x, c_y, c_z)$, radius $= \sqrt{k + |c|^2}$.
- **Gauss-Newton refinement** (when ≥6 vertices): After the algebraic fit, iteratively minimize the geometric distance $f_i = \|p_i - c\| - r$. The Jacobian is $J_i = [-(p_i - c)/\|p_i - c\|, -1]$. Solve the normal equations $J^T J \delta = -J^T f$ at each iteration (4×4 system via Gaussian elimination). Up to 20 iterations, converges when $\|\delta\| < 10^{-14}$. This refines the algebraic fit to minimize the true vertex-to-surface distance rather than the algebraic residual.

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

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. ~~Key challenge: partial spherical patches are hard to distinguish from cylinders or even planes if the patch is small relative to the radius. Consider requiring a minimum angular extent or using curvature analysis to disambiguate.~~ Addressed: solid-angle coverage validation prevents 1D strip growth (cylinder fillets), surface-tolerance validation in BFS prevents growth along non-spherical surfaces, and max_sphere_radius prevents degenerate large-radius fits on flat faces. Also consider toroidal surfaces (donuts, fillets)—they're common in CAD and combine characteristics of cylinders and spheres. Note: torus faces locally fit spheres, producing sphere hypotheses on torus surfaces (e.g., pipe_elbow). Stage 2.6 will need torus hypothesis support to resolve this.*

---

## 2.4 Deduce conical hypotheses

Fit conical hypotheses to connected sets of faces that lie on a common cone. A cone is a surface of revolution (SOR) with a linear profile — it generalizes the cylinder (which has a constant-radius profile). Conical surfaces commonly appear as chamfers, tapered bores, and funnel shapes in CAD models.

A conical hypothesis consists of:
- `apex: [f64; 3]` — Apex (tip) of the cone.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis (from apex outward toward the open end).
- `half_angle: f64` — Half-angle of the cone (in radians). The angle between the axis and any generator line on the surface.
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a cone) or toward it (concave, like the inside of a conical hole).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this cone.
- `error_max: f64` — Maximum absolute distance from any vertex to the cone surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

**Cone parameterization:**
A cone with apex $A$, axis direction $\hat{a}$, and half-angle $\theta$ satisfies: for any point $p$ on the surface, the angle between $(p - A)$ and $\hat{a}$ equals $\theta$. Equivalently, in (h,r) coordinates where $h = (p - A) \cdot \hat{a}$ (signed axial distance from apex) and $r = |(p - A) - h\hat{a}|$ (radial distance from axis), the surface satisfies $r = h \tan\theta$.

The signed distance from a point $p$ to the cone surface (for LM residual) is:
$$d = r_i \cos\theta - h_i \sin\theta$$
where $h_i = (p_i - A) \cdot \hat{a}$ and $r_i = |(p_i - A) - h_i \hat{a}|$.

**Initialization via (h,r) profile fitting:**
This exploits the key insight from the unified SOR framework (see Algorithm Reference below): all surfaces of revolution reduce to a 2D profile fitting problem once the axis is estimated.

1. **Axis estimation** from normal covariance (same as cylinder): compute weighted covariance $M = \sum w_i n_i n_i^T$ where $n_i$ are face normals and $w_i$ are face areas. The eigenvector with smallest eigenvalue is the axis direction estimate.
2. **Reduce to (h,r) coordinates**: $h_i = (v_i - \text{centroid}) \cdot \hat{a}$, $r_i = |(v_i - \text{centroid}) - h_i \hat{a}|$.
3. **Fit linear profile** $r = m \cdot h + b$ via weighted linear regression (weights = vertex area).
4. **Extract cone parameters**: half-angle $\theta = \arctan(m)$, apex at distance $-b/m$ from centroid along axis.
5. **Reject if cylindrical**: if $|m| < 0.01$ (nearly constant radius), this is a cylinder, not a cone. Also reject if $\theta > 85°$ (nearly planar).

**Alternative initialization via quadric fitting:**
A cone is a quadric surface: $x^T A x + b^T x + c = 0$ with specific eigenvalue structure. This approach is useful as a cross-check or when the axis is poorly estimated.

1. Build the 10-coefficient quadric feature matrix: for each vertex $p_i = (x,y,z)$, construct row $[x^2, y^2, z^2, xy, xz, yz, x, y, z, 1]$.
2. SVD of the feature matrix. The right singular vector with smallest singular value gives the quadric coefficients.
3. Convert to matrix form: $Q(x) = x^T \mathbf{A} x + b^T x + c$.
4. Extract cone parameters:
   - Apex: $A = -\frac{1}{2} \mathbf{A}^{-1} b$ (gradient of Q is zero at apex).
   - Axis: eigenvector of $\mathbf{A}$ corresponding to the distinct eigenvalue (for a cone, two eigenvalues are equal and one is different).
   - Half-angle: from the eigenvalue ratio.
5. Validate: eigenvalue pattern must match a cone (two eigenvalues approximately equal, one distinct, with specific sign relationships). If the pattern matches a cylinder (one eigenvalue ≈ 0), sphere (all equal), or ellipsoid, skip.

**Apex-vertex seeding and BFS region growing:**
- Identify apex candidates: vertices with ≥6 incident faces that have non-coplanar normals (smallest eigenvalue / trace < 0.3 from normal covariance). Skip faces already committed to cylindrical hypotheses.
- For each apex candidate, use ALL incident faces as the seed set. Fit the cone from the seed, then early-reject if half-angle is outside [2°, 85°]. Verify all seed vertices within vertex_tol of the fitted cone (with multi-axis LM rescue if the initial fit doesn't reach vertex_tol).
- Expand via BFS: for each candidate neighbor face, check vertex-to-cone distance < vertex_tol. If within REFIT_SKIP_MULTIPLIER × vertex_tol but not vertex_tol, re-fit and accept only if the refit satisfies vertex_tol (all vertices) and surface_tol (all centroids). Skip faces already committed to conical or cylindrical hypotheses.
- Centroid validation: face centroid within surface_tolerance of the cone.
- Angular tolerance: dihedral angle between adjacent faces ≤ 2× angular_tolerance (relaxed because tessellated cones have wider inter-strip angles than cylinders).
- Angular coverage validation (same algorithm as cylinder).
- Minimum 6-face requirement (a cone has 5 independent DOF but 6 faces provides margin).

**Distinguishing cones from cylinders:**
After fitting both a cylinder and a cone to the same face set, compare RMS residuals. If the cone fit is not significantly better than the cylinder fit (< 20% RMS improvement), prefer the simpler cylinder interpretation. This prevents noisy flat-profile data from being misclassified as a very gentle cone.

**LM refinement:**
Optimize $[A_x, A_y, A_z, \alpha, \beta, \theta]$ (6 parameters: apex position 3, axis tilt angles 2, half-angle 1).
Residual: $d_i = r_i \cos\theta - h_i \sin\theta$.

- With the `--compare` flag:
    - Same centroid-projection validation as cylinder/sphere.

*Comment: The (h,r) profile approach is elegant because it handles the cylinder/cone ambiguity naturally: a horizontal line in (h,r) is a cylinder, a sloped line is a cone, and the distinction is just the slope parameter. This is described in the SOR reduction technique. The quadric approach provides an independent validation path.*

---

## 2.5 Deduce toroidal hypotheses

Fit toroidal hypotheses to connected sets of faces that lie on a common torus. Toroidal surfaces commonly appear as fillets (rounded edges), blends, and pipe elbows in CAD models. A torus is a surface of revolution with a circular profile offset from the axis.

A toroidal hypothesis consists of:
- `center: [f64; 3]` — Center of the torus (on the axis, in the plane of the major circle).
- `axis_direction: [f64; 3]` — Unit direction vector along the torus axis.
- `major_radius: f64` — Distance from center to the tube centerline ($R$).
- `minor_radius: f64` — Radius of the tube cross-section ($r$). For a fillet of radius $r$, the minor radius is $r$.
- `convex: bool` — Whether the torus patch is on the outer (convex) or inner (concave) surface of the tube.
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this torus.
- `error_max: f64` — Maximum absolute distance from any vertex to the torus surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

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
- LM patience 500, minimum 7 data points. Numerical Jacobian via central differences (eps=1e-8).

**Alternative fitting: SOR axis method (`fit_torus_sor`).**
When the tube center method fails (e.g., too few tube centers above threshold or quality gate failure), a surface-of-revolution axis estimation is used as fallback. Estimates torus axis from the smallest eigenvector of normal covariance, projects vertices to (h,r) coordinates, fits a shifted circle to deduce R and r. Same LM refinement follows. `fit_torus` tries **both** methods (tube center and SOR) and keeps whichever achieves the lowest maximum vertex-to-torus distance.

**Deterministic iteration.**
All vertex set iterations (HashSet) are sorted before use to ensure deterministic results regardless of hash ordering. Without this, the fitting results can vary between runs due to floating-point non-associativity.

**Seeding strategy:**
Face-based seeding: for each uncommitted candidate face (normal, non-hypothesized), build a seed consisting of that face plus all faces that share at least one vertex with it (the **1-ring vertex neighborhood**). Only uncommitted faces are included in the seed. The seed is sorted by face index for deterministic results. Seeds with fewer than `MIN_TORUS_FACES` (7) faces or fewer than 9 unique vertices are skipped. BFS skips faces that already have a toroidal hypothesis, so later seeds naturally avoid regions already identified.

**Multi-axis LM rescue:**
If the initial geometric fit (tube center or SOR) does not meet vertex_tolerance, a multi-axis Levenberg–Marquardt rescue is attempted. Five axis candidates are tried: the geometric estimate, the smallest eigenvector of the face-normal covariance matrix, and the three coordinate axes [1,0,0], [0,1,0], [0,0,1]. For each axis candidate:
1. Project vertices to the plane perpendicular to the axis and fit a 2D algebraic circle (Kasa method, 3×3 linear system) to estimate center position and major radius $R$.
2. Estimate minor radius $r$ from the residuals of the 2D circle fit.
3. Run LM refinement from this initial estimate.
4. Track the best result across all axis candidates. Accept if the best achieves vertex_tolerance.

**BFS region growing:**
- Vertex-to-torus distance: $|\sqrt{(\text{radial} - R)^2 + \text{axial}^2} - r| < \text{vertex\_tolerance}$.
- Centroid validation: face centroid within surface_tolerance of torus.
- Convexity check: normal direction consistent with hypothesis.
- Angular tolerance: dihedral angle ≤ angular_tolerance.
- Re-fit on marginal faces: if a face's vertices are within `REFIT_SKIP_MULTIPLIER × vertex_tolerance`, attempt a refit with the face included; accept iff the refit satisfies both vertex_tolerance (all vertices) and surface_tolerance (all centroids).

**Post-growth validation:**
- Minimum 7-face requirement (torus has 7 DOF, LM needs ≥7 points).
- Angular coverage around the minor circle: faces should span multiple angular positions around the tube cross-section.
- Final vertex_tolerance re-validation after the last re-fit.
- Face centroid within surface_tolerance.

*Comment: The medial axis / tube center method is far more robust for fillet-sized torus patches than direct 7-parameter nonlinear fitting. By reducing the problem to (1) estimate one scalar r, (2) fit a 3D circle to the tube centers, it separates the easy part (circle fitting is linear algebra) from the hard part (estimating the minor radius). The spine classification step (Step 3) also provides a natural way to detect degenerate cases where the "torus" is actually a cylinder or sphere — this happens when the fillet radius approaches the edge length or when the surface is nearly spherical. The face-based seeding strategy (one seed per face) is simpler and more robust than vertex-based seeding with post-merge: every toroidal face gets a chance to anchor a hypothesis, and BFS naturally deduplicates by skipping already-committed faces.*

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

The current implementation uses separate fitting algorithms for each surface type (stages 2.2–2.5), which is appropriate for incremental development. A future refactoring pass (see Implementation Phases in DEVELOPMENT_PLAN.md) will unify these into the SOR framework to improve accuracy on ambiguous cases and simplify the codebase.

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

This is used in stage 2.5 (torus fitting): estimate $r$ from normal-line intersections between adjacent faces, compute tube centers $k_i$, fit a 3D circle to get the major circle, and extract the torus parameters.

### Quadric Fitting

Any quadric surface (plane, sphere, cylinder, cone, ellipsoid, hyperboloid, paraboloid) can be written as:
$$Q(x, y, z) = Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0$$

Fitting: construct the $N \times 10$ feature matrix and find its null space via SVD (smallest singular vector). The matrix form $x^T \mathbf{A} x + b^T x + c = 0$ encodes the surface type in the eigenvalues of $\mathbf{A}$:
- **Sphere**: $\mathbf{A} = \lambda I$ (three equal eigenvalues).
- **Cylinder**: one eigenvalue is 0 (axis direction), other two equal.
- **Cone**: three nonzero eigenvalues with $\det(\mathbf{A}) \approx 0$ (passes through apex).
- **Plane**: $\mathbf{A} = 0$ (all eigenvalues zero).

Quadric fitting is fast (SVD of an $N \times 10$ matrix, $O(N)$) and provides a good initialization for specialized fitters. It's used as an alternative initialization for cone fitting (stage 2.4).
