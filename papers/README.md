# Research Papers: Fitting Analytical Surfaces to CAD-Generated Tessellations

Literature survey on fitting analytical curved surfaces (planes, cylinders, spheres, cones, tori) to tessellated meshes where vertices lie on the approximated surfaces. Organized by research group and approach.

## Downloaded Papers

| File | Paper |
|------|-------|
| `Schnabel_2007_Efficient_RANSAC.pdf` | Schnabel, Wahl, Klein 2007 — Efficient RANSAC |
| `Benko_2001_RE_BRep_Algorithms.pdf` | Benkő, Martin, Várady 2001 — RE B-Rep Algorithms |
| `Benko_2002_Constrained_Fitting.pdf` | Benkő, Kós, Várady, Martin 2002 — Constrained Fitting |
| `Fit4CAD_2022_Benchmark.pdf` | Romanengo et al. 2022 — Fit4CAD Benchmark |

**Not downloaded** (behind paywalls, no open-access version found):
- Várady, Martin, Cox 1997 — behind Elsevier paywall (DOI: 10.1016/S0010-4485(96)00054-1)
- Benkő, Várady 2004 — behind Elsevier paywall (DOI: 10.1016/S0010-4485(03)00160-2)
- GlobFit (Li et al. 2011) — behind ACM paywall (DOI: 10.1145/1964921.1964947)
- SHREC 2022 — behind paywall

## Relevance to Brepper

Most reverse engineering literature addresses **noisy scanned point clouds** (laser scanners, structured light), where measurement noise is the primary challenge. Brepper's problem is different: vertices lie **exactly** (within floating-point epsilon) on the original CAD surfaces. The challenge is not noise rejection but **structural outliers** — faces at surface boundaries that straddle two surfaces, and the combinatorial problem of finding which faces belong to which surface.

Key algorithmic ideas borrowed from this literature:
- **Region growing** from seed faces (Benkő/Várady) — directly used in brepper's BFS expansion
- **Normal-based axis estimation** for cylinders (eigenvector of normal covariance matrix) — used in brepper's cylinder fitting
- **Gaussian image** analysis for cylinder axis candidates (Schnabel) — related to brepper's cross-product seeding
- **Multi-seed evaluation** to avoid committing to suboptimal seeds — brepper's enhancement over the basic region growing approach

## Foundational Work: Várady/Benkő Group (Budapest)

### Várady, Martin, Cox — "Reverse engineering of geometric models" (1997)
- **Journal**: Computer-Aided Design, Vol. 29, No. 4, pp. 255–268
- **Cited**: ~2200 times
- **Summary**: Foundational survey of the entire reverse engineering pipeline: data acquisition → preprocessing → segmentation → surface fitting → CAD model construction. Defines the taxonomy of approaches (region growing vs. edge detection vs. hybrid) that subsequent work follows.
- **DOI**: 10.1016/S0010-4485(96)00054-1

### Benkő, Martin, Várady — "Algorithms for reverse engineering boundary representation models" (2001)
- **Journal**: Computer-Aided Design, Vol. 33, No. 11, pp. 839–851
- **Cited**: ~340 times
- **Summary**: Full pipeline from tessellated mesh to B-Rep: segmentation via region growing, fitting of planes/cylinders/spheres/cones/tori, constrained refitting, topology construction. Closest prior art to brepper's overall approach. Uses curvature estimation at vertices to classify surface type before fitting.
- **DOI**: 10.1016/S0010-4485(01)00012-0

### Benkő, Várady — "Segmentation methods for smooth point regions of conventional engineering objects" (2004)
- **Journal**: Computer-Aided Design, Vol. 36, No. 6, pp. 511–523
- **Cited**: ~155 times
- **Summary**: Compares three segmentation approaches: (1) region growing from seed triangles, (2) direct decomposition using curvature sign changes, (3) hybrid. Concludes region growing with good seed selection outperforms pure boundary detection. Directly relevant to brepper's choice of BFS region growing.
- **DOI**: 10.1016/S0010-4485(03)00160-2

### Benkő, Kós, Várady, Martin — "Constrained fitting in reverse engineering" (2002)
- **Journal**: Computer Aided Geometric Design, Vol. 19, No. 3, pp. 173–205
- **Cited**: ~280 times
- **Summary**: After initial fitting, applies geometric constraints (parallelism, perpendicularity, coaxiality, tangency) to improve the model. E.g., two cylinders that are "almost coaxial" are constrained to share an axis. Relevant to brepper's future surface refitting stage (2.7).
- **DOI**: 10.1016/S0167-8396(01)00085-1

## RANSAC Approaches

### Schnabel, Wahl, Klein — "Efficient RANSAC for Point-Cloud Shape Detection" (2007)
- **Journal**: Computer Graphics Forum, Vol. 26, No. 2, pp. 214–226
- **Cited**: ~2900 times
- **Summary**: Adapts RANSAC to detect multiple primitive shapes (planes, cylinders, spheres, cones, tori) simultaneously in unorganized point clouds. Key innovations: (1) uses **Gaussian image** (normal sphere) to efficiently generate cylinder/cone axis candidates, (2) bitmap-based connected component validation, (3) processes shapes largest-first. The multi-seed evaluation approach in brepper is conceptually related — try many seeds, keep the best — though brepper uses exhaustive enumeration of neighbor pairs rather than random sampling.
- **DOI**: 10.1111/j.1467-8659.2007.01016.x
- **Note**: Open-source implementation available (C++); GlobFit (Li et al. 2011) extends this with global constraint optimization.

## Benchmarks and Datasets

### Kaiser, Dalstein, Oesau, Lafarge — "Fit4CAD: A point cloud benchmark for fitting simple primitives" (2022)
- **Conference**: ECCV 2022 Workshop
- **Cited**: ~35 times
- **Summary**: Standardized benchmark for evaluating primitive fitting on CAD-like point clouds. Includes ground truth segmentation and fitted parameters for planes, cylinders, spheres, cones. Useful for comparing brepper's accuracy against other methods.

### SHREC 2022 — "Fitting and recognition of geometric primitives in segmented 3D point clouds" (2022)
- **Cited**: ~16 times
- **Summary**: Shape retrieval contest track specifically for primitive fitting. Multiple teams competed with different algorithms, providing a comparison of state-of-the-art approaches.

## Additional References

### Li, Wu, Sharf, Cohen-Or, Chen — "GlobFit: Consistently Fitting Primitives by Discovering Global Relations" (2011)
- **Conference**: ACM SIGGRAPH
- **Summary**: Extends Schnabel's RANSAC with global relation detection (parallelism, orthogonality, coplanarity, equal radius/spacing) and simultaneous constraint optimization. Relevant to brepper's future constrained refitting.
- **DOI**: 10.1145/1964921.1964947

### Attene, Falcidieno, Spagnuolo — "Hierarchical mesh segmentation based on fitting primitives" (2006)
- **Journal**: The Visual Computer, Vol. 22, No. 3, pp. 181–193
- **Summary**: Hierarchical approach: start with each triangle as a region, merge adjacent regions if they fit the same primitive. Bottom-up alternative to brepper's top-down seed-and-grow approach.
