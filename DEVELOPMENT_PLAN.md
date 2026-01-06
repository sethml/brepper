# STL to STEP Conversion Utility - Development Plan

## Project Overview

**brepper** (B-Rep from Mesh) - A C++ utility to convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Dependencies

| Library | Purpose | Version |
|---------|---------|---------|
| PCL (Point Cloud Library) | Mesh I/O, point cloud processing, RANSAC segmentation | ≥1.12 |
| OpenCASCADE (OCCT) | B-Rep modeling, surface fitting, STEP export | ≥7.6 |
| Eigen | Linear algebra (bundled with PCL) | ≥3.4 |
| CLI11 or cxxopts | Command-line argument parsing | latest |

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              Pipeline Stages                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐  │
│  │  Stage 1 │ → │  Stage 2 │ → │  Stage 3 │ → │  Stage 4 │ → │  Stage 5 │  │
│  │  Mesh    │   │  Point   │   │  Surface │   │  Boundary│   │  B-Rep   │  │
│  │  Input   │   │  Cloud   │   │  Fitting │   │  Detect  │   │  Output  │  │
│  └──────────┘   └──────────┘   └──────────┘   └──────────┘   └──────────┘  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Stage 1: Mesh Input & Preprocessing

### 1.1 Read STL File
- **Library**: `pcl::io::loadPolygonFileSTL()` or `pcl::io::loadPolygonFile()`
- **Input**: Binary or ASCII STL file
- **Output**: `pcl::PolygonMesh`

### 1.2 Compute Triangle Normals
- Compute face normals for each triangle: `n = normalize((v1-v0) × (v2-v0))`
- Store face normals in a parallel data structure
- Optionally compute vertex normals by averaging adjacent face normals

### 1.3 Uniform Triangle Sampling
- Sample additional points within triangles for denser point clouds
- Use barycentric coordinate sampling: `P = α·v0 + β·v1 + γ·v2` where `α+β+γ=1`
- Sampling density should be configurable (points per unit area or per triangle)

**Parameters:**
- `--sample-density <float>` - Points per unit area (default: auto-compute based on mesh resolution)
- `--min-samples-per-triangle <int>` - Minimum samples per triangle (default: 1)

### 1.4 Convert to Point Cloud with Normals
- **Output**: `pcl::PointCloud<pcl::PointNormal>`
- Interpolate normals at sampled points from vertex normals or use face normal

---

## Stage 2: Surface Segmentation (Iterative RANSAC)

### 2.1 RANSAC Surface Fitting Loop

```cpp
while (remaining_points > threshold) {
    // 1. Run SAC segmentation for each model type
    // 2. Select best fit (highest inlier count meeting quality threshold)
    // 3. Extract inliers, add to surface list
    // 4. Remove inliers from point cloud
}
```

### 2.2 Supported Surface Types (Priority Order)

| Surface Type | PCL Model | OCCT Surface |
|--------------|-----------|--------------|
| Plane | `SACMODEL_PLANE` | `Geom_Plane` |
| Cylinder | `SACMODEL_CYLINDER` | `Geom_CylindricalSurface` |
| Sphere | `SACMODEL_SPHERE` | `Geom_SphericalSurface` |
| Cone | `SACMODEL_CONE` | `Geom_ConicalSurface` |
| Torus | Custom implementation | `Geom_ToroidalSurface` |

### 2.3 PCL Segmentation Setup

```cpp
pcl::SACSegmentationFromNormals<PointT, NormalT> seg;
seg.setOptimizeCoefficients(true);
seg.setMethodType(pcl::SAC_RANSAC);
seg.setMaxIterations(max_iterations);
seg.setDistanceThreshold(distance_threshold);
seg.setNormalDistanceWeight(normal_weight);
```

### 2.4 Extract Inliers

```cpp
pcl::ExtractIndices<PointT> extract;
extract.setInputCloud(cloud);
extract.setIndices(inliers);
extract.setNegative(false);  // Get inliers
extract.filter(*surface_cloud);
extract.setNegative(true);   // Get remaining
extract.filter(*remaining_cloud);
```

**Parameters:**
- `--ransac-distance <float>` - RANSAC distance threshold (default: 0.01)
- `--ransac-iterations <int>` - Max RANSAC iterations (default: 1000)
- `--normal-weight <float>` - Normal consistency weight (default: 0.1)
- `--min-inliers <int>` - Minimum inliers for valid surface (default: 100)
- `--min-inlier-ratio <float>` - Minimum ratio of inliers (default: 0.01)

---

## Stage 3: Clustering & NURBS Fitting

### 3.1 Cluster Similar Primitive Segments

Merge segments of the same type with similar parameters:
- Planes with similar normals and close distances
- Cylinders with similar axes and radii
- etc.

**Parameters:**
- `--plane-angle-threshold <float>` - Max angle between plane normals to merge (degrees, default: 5.0)
- `--plane-distance-threshold <float>` - Max distance between planes to merge (default: 0.01)
- `--cylinder-radius-threshold <float>` - Max radius difference to merge (default: 0.01)

### 3.2 Euclidean Clustering for Remaining Points

```cpp
pcl::EuclideanClusterExtraction<PointT> ec;
ec.setClusterTolerance(cluster_tolerance);
ec.setMinClusterSize(min_cluster_size);
ec.setMaxClusterSize(max_cluster_size);
ec.setInputCloud(remaining_cloud);
ec.extract(cluster_indices);
```

**Parameters:**
- `--cluster-tolerance <float>` - Clustering distance (default: 0.02)
- `--min-cluster-size <int>` - Minimum cluster points (default: 50)
- `--max-cluster-size <int>` - Maximum cluster points (default: 1000000)

### 3.3 Fit NURBS Surfaces to Clusters

For remaining organic/freeform regions:
- Use OCCT's `GeomAPI_PointsToBSplineSurface`
- Or PCL's B-spline fitting: `pcl::on_nurbs::FittingSurface`

**Parameters:**
- `--nurbs-degree <int>` - NURBS surface degree (default: 3)
- `--nurbs-control-points <int>` - Control points per direction (default: 10)
- `--nurbs-fitting-tolerance <float>` - Fitting tolerance (default: 0.001)

---

## Stage 4: Mesh Segmentation & Boundary Detection

### 4.1 Assign Triangles to Surfaces

For each mesh triangle:
1. Compute centroid and normal
2. Find closest fitted surface
3. Check distance and normal consistency
4. Assign triangle to surface ID

```cpp
struct TriangleAssignment {
    int triangle_id;
    int surface_id;
    double distance;
    double normal_deviation;
};
```

**Parameters:**
- `--assignment-distance <float>` - Max distance for triangle assignment (default: 0.02)
- `--assignment-angle <float>` - Max normal deviation (degrees, default: 15.0)

### 4.2 Detect Boundary Edges

```cpp
for each edge (v0, v1) in mesh:
    tri1 = adjacent_triangle_1(edge)
    tri2 = adjacent_triangle_2(edge)
    if surface_id[tri1] != surface_id[tri2]:
        boundary_edges.add(edge, surface_id[tri1], surface_id[tri2])
```

### 4.3 Group Edges into Boundary Curves

1. Build edge adjacency graph
2. Extract connected chains of boundary edges
3. Order vertices along each chain

```cpp
struct BoundaryCurve {
    std::vector<Eigen::Vector3d> points;
    int surface_id_left;
    int surface_id_right;
};
```

---

## Stage 5: Curve Fitting & B-Rep Construction

### 5.1 Fit Curves to Boundary Chains

| Surface Pair | Curve Type | OCCT Class |
|--------------|------------|------------|
| Plane-Plane | Line | `Geom_Line` |
| Plane-Cylinder | Line or Ellipse | `Geom_Line` / `Geom_Ellipse` |
| Cylinder-Cylinder | Line or Ellipse | varies |
| Plane-Sphere | Circle | `Geom_Circle` |
| Any-NURBS | B-Spline | `Geom_BSplineCurve` |
| Other | B-Spline | `Geom_BSplineCurve` |

**Implementation:**
```cpp
// For analytic curves
Handle(Geom_Curve) fitBoundaryCurve(
    const BoundaryCurve& boundary,
    SurfaceType type1, 
    SurfaceType type2
);
```

**Parameters:**
- `--curve-fitting-tolerance <float>` - Curve fitting tolerance (default: 0.001)
- `--prefer-analytic-curves <bool>` - Prefer lines/circles over splines (default: true)

### 5.2 Create OCCT Surfaces

Convert fitted surfaces to OCCT:

```cpp
Handle(Geom_Surface) toOCCT(const FittedSurface& surface) {
    switch (surface.type) {
        case PLANE:
            return new Geom_Plane(gp_Pln(origin, normal));
        case CYLINDER:
            return new Geom_CylindricalSurface(gp_Ax3(origin, axis), radius);
        // ... etc
    }
}
```

### 5.3 Trim Surfaces with Curves

1. Create `TopoDS_Edge` from each boundary curve
2. Create `TopoDS_Wire` from connected edges
3. Create `TopoDS_Face` from surface + wire boundary

```cpp
BRepBuilderAPI_MakeFace faceMaker(surface, wire, true);
TopoDS_Face face = faceMaker.Face();
```

### 5.4 Heal and Sew Faces into Shell

```cpp
BRepBuilderAPI_Sewing sewing(tolerance);
for (auto& face : faces) {
    sewing.Add(face);
}
sewing.Perform();
TopoDS_Shape shell = sewing.SewedShape();
```

**Parameters:**
- `--sewing-tolerance <float>` - Sewing tolerance (default: 0.001)

### 5.5 Create Solid from Shell

```cpp
BRepBuilderAPI_MakeSolid solidMaker;
solidMaker.Add(TopoDS::Shell(shell));
TopoDS_Solid solid = solidMaker.Solid();

// Fix orientation
BRepLib::OrientClosedSolid(solid);
```

### 5.6 Shape Healing

```cpp
ShapeFix_Shape fixer(solid);
fixer.SetPrecision(tolerance);
fixer.Perform();
TopoDS_Shape fixed = fixer.Shape();
```

**Parameters:**
- `--healing-tolerance <float>` - Shape healing tolerance (default: 0.001)

---

## Stage 6: STEP Export

```cpp
STEPControl_Writer writer;
writer.Transfer(solid, STEPControl_AsIs);
IFSelect_ReturnStatus status = writer.Write(output_path);
```

**Parameters:**
- `--step-schema <string>` - STEP schema: AP203, AP214, AP242 (default: AP214)

---

## Command Line Interface

```
brepper - Convert STL mesh to STEP with fitted surfaces

USAGE:
    brepper [OPTIONS] <input.stl> -o <output.step>

REQUIRED:
    <input.stl>              Input STL file (binary or ASCII)
    -o, --output <file>      Output STEP file

GENERAL OPTIONS:
    -v, --verbose            Enable verbose output
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --threads <N>            Number of threads (default: auto)

MESH PREPROCESSING:
    --sample-density <F>     Points per unit area (default: auto)
    --min-samples <N>        Min samples per triangle (default: 1)

RANSAC SEGMENTATION:
    --ransac-distance <F>    Distance threshold (default: 0.01)
    --ransac-iterations <N>  Max iterations (default: 1000)
    --normal-weight <F>      Normal weight 0-1 (default: 0.1)
    --min-inliers <N>        Min points per surface (default: 100)
    --min-inlier-ratio <F>   Min ratio of cloud (default: 0.01)

SURFACE TYPES:
    --fit-planes             Fit planes (default: on)
    --fit-cylinders          Fit cylinders (default: on)
    --fit-spheres            Fit spheres (default: on)
    --fit-cones              Fit cones (default: on)
    --fit-tori               Fit tori (default: off)
    --no-<type>              Disable specific surface type

CLUSTERING:
    --plane-merge-angle <F>  Plane merge angle threshold (deg, default: 5.0)
    --plane-merge-dist <F>   Plane merge distance (default: 0.01)
    --cluster-tolerance <F>  Euclidean cluster tolerance (default: 0.02)
    --min-cluster-size <N>   Minimum cluster size (default: 50)

NURBS FITTING:
    --nurbs-degree <N>       B-spline degree (default: 3)
    --nurbs-refinement <N>   Control point density (default: 10)
    --nurbs-tolerance <F>    Fitting tolerance (default: 0.001)

TRIANGLE ASSIGNMENT:
    --assign-distance <F>    Max assignment distance (default: 0.02)
    --assign-angle <F>       Max normal deviation (deg, default: 15.0)

CURVE FITTING:
    --curve-tolerance <F>    Curve fitting tolerance (default: 0.001)
    --prefer-analytic        Prefer analytic curves (default: true)

B-REP CONSTRUCTION:
    --sewing-tolerance <F>   Face sewing tolerance (default: 0.001)
    --healing-tolerance <F>  Shape healing tolerance (default: 0.001)
    --step-schema <S>        STEP schema: AP203|AP214|AP242 (default: AP214)

DEBUG OUTPUT:
    --save-point-cloud <F>   Save sampled point cloud (PCD/PLY)
    --save-segmentation <F>  Save segmented mesh (PLY with colors)
    --save-boundaries <F>    Save boundary curves (PLY)

PRESETS:
    --preset <name>          Use parameter preset:
                             - tight: High precision CAD export
                             - loose: Low-quality mesh repair
                             - default: Balanced settings
```

---

## File Structure

```
brepper/
├── CMakeLists.txt
├── README.md
├── DEVELOPMENT_PLAN.md
├── src/
│   ├── main.cpp                    # Entry point, CLI parsing
│   ├── brepper.hpp                 # Main pipeline orchestration
│   ├── brepper.cpp
│   ├── mesh/
│   │   ├── stl_reader.hpp          # STL file reading
│   │   ├── stl_reader.cpp
│   │   ├── mesh_sampling.hpp       # Triangle sampling
│   │   ├── mesh_sampling.cpp
│   │   ├── normal_computation.hpp  # Normal computation
│   │   └── normal_computation.cpp
│   ├── segmentation/
│   │   ├── ransac_segmenter.hpp    # RANSAC surface fitting
│   │   ├── ransac_segmenter.cpp
│   │   ├── surface_clustering.hpp  # Segment merging
│   │   ├── surface_clustering.cpp
│   │   ├── nurbs_fitter.hpp        # NURBS surface fitting
│   │   └── nurbs_fitter.cpp
│   ├── boundary/
│   │   ├── triangle_assignment.hpp # Triangle-to-surface mapping
│   │   ├── triangle_assignment.cpp
│   │   ├── edge_detection.hpp      # Boundary edge detection
│   │   ├── edge_detection.cpp
│   │   ├── curve_extraction.hpp    # Edge chain grouping
│   │   ├── curve_extraction.cpp
│   │   ├── curve_fitting.hpp       # Analytic curve fitting
│   │   └── curve_fitting.cpp
│   ├── brep/
│   │   ├── surface_converter.hpp   # PCL → OCCT surface
│   │   ├── surface_converter.cpp
│   │   ├── face_builder.hpp        # Trimmed face construction
│   │   ├── face_builder.cpp
│   │   ├── shell_builder.hpp       # Sewing & shell construction
│   │   ├── shell_builder.cpp
│   │   ├── solid_builder.hpp       # Solid creation & healing
│   │   └── solid_builder.cpp
│   ├── io/
│   │   ├── step_writer.hpp         # STEP export
│   │   └── step_writer.cpp
│   └── common/
│       ├── types.hpp               # Common data structures
│       ├── config.hpp              # Configuration parameters
│       └── logging.hpp             # Logging utilities
├── tests/
│   ├── test_stl_reader.cpp
│   ├── test_ransac.cpp
│   ├── test_boundary.cpp
│   └── test_data/
│       ├── cube.stl
│       ├── cylinder.stl
│       └── complex.stl
└── examples/
    ├── simple_cube.stl
    └── README.md
```

---

## Implementation Phases

### Phase 1: Foundation (Week 1-2)
- [x] Project setup (CMake, dependencies)
- [x] STL file reading
- [x] Point cloud generation with normals
- [x] Basic CLI framework

### Phase 2: Surface Segmentation (Week 3-4)
- [x] RANSAC plane fitting
- [x] RANSAC cylinder fitting
- [x] RANSAC sphere/cone fitting
- [x] Iterative extraction loop
- [ ] Segment clustering

### Phase 3: NURBS & Assignment (Week 5-6)
- [ ] Euclidean clustering
- [ ] NURBS surface fitting
- [x] Triangle-to-surface assignment
- [ ] Segmented mesh visualization

### Phase 4: Boundary Detection (Week 7-8)
- [ ] Boundary edge detection
- [ ] Edge chain extraction
- [ ] Analytic curve fitting
- [ ] B-spline curve fitting

### Phase 5: B-Rep Construction (Week 9-10)
- [ ] OCCT surface creation
- [ ] Trimmed face construction
- [ ] Face sewing
- [ ] Solid creation & healing

### Phase 6: Export & Polish (Week 11-12)
- [ ] STEP export
- [ ] Error handling & validation
- [x] Performance optimization (OpenMP parallelization)
- [x] Testing & documentation (Catch2 unit tests)

---

## Key Challenges & Mitigations

| Challenge | Mitigation Strategy |
|-----------|---------------------|
| Noisy STL meshes | Robust RANSAC, adjustable thresholds |
| Degenerate triangles | Pre-filter invalid geometry |
| Ambiguous surface types | Score-based selection, user hints |
| Boundary curve discontinuities | Smoothing, tolerance handling |
| Topology errors in B-Rep | OCCT ShapeFix, iterative healing |
| Performance with large meshes | Spatial indexing (KD-tree), parallelization |

---

## Testing Strategy

1. **Unit Tests**: Individual components (readers, fitters, converters)
2. **Integration Tests**: Full pipeline on known geometries
3. **Regression Tests**: Compare output STEP to reference
4. **Validation**: Round-trip test (STEP → mesh → STEP)

### Test Cases

| Test Case | Input | Expected Output |
|-----------|-------|-----------------|
| Cube | 12 triangles | 6 planes |
| Cylinder | Tessellated cylinder | 1 cylinder + 2 planes |
| Sphere | Tessellated sphere | 1 sphere |
| Fillet | Blended edge | Planes + fillet surface |
| Complex part | Real CAD export | Matching topology |

---

## Performance Targets

| Mesh Size | Target Time |
|-----------|-------------|
| < 10K triangles | < 5 seconds |
| 10K - 100K triangles | < 30 seconds |
| 100K - 1M triangles | < 5 minutes |

---

## Future Enhancements

- [ ] Support for OBJ, PLY input formats
- [ ] IGES export option
- [ ] GUI for interactive parameter tuning
- [ ] Machine learning for surface type classification
- [ ] Hole detection and filling
- [ ] Feature recognition (holes, pockets, bosses)
