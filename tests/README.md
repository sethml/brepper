# Test Data for brepper

This directory contains test STL files organized by source:

## Directory Structure

- `manual/` - Hand-written ASCII STL files for basic testing
- `generated/` - STL files exported from CAD programs (Fusion 360, Onshape)

## Manual Test Files

### cube.stl
A unit cube (0,0,0) to (1,1,1) with 12 triangles. Used for basic import and sampling tests.

## Recommended Test Models to Export

When exporting from Fusion 360 or Onshape, use the following models to test different surface types and mesh densities:

### Basic Geometric Primitives
1. **sphere_25mm.stl** - 25mm diameter sphere
   - Tests: spherical surface fitting, uniform curvature
   - Export at: Fine mesh (many small triangles)

2. **cylinder_10x30mm.stl** - 10mm diameter, 30mm height cylinder
   - Tests: cylindrical surface fitting, flat end caps
   - Export at: Normal mesh

3. **cone_15x20mm.stl** - 15mm base diameter, 20mm height cone
   - Tests: conical surface fitting, point singularity
   - Export at: Normal mesh

### Multi-Surface Models
4. **rounded_cube_10mm_r2.stl** - S cube with 2mm edge fillets
   - Tests: plane + cylindrical fillet detection, surface transitions
   - Export at: Fine mesh to capture fillets

5. **pipe_elbow_10mm.stl** - 90° pipe elbow, 10mm inner diameter
   - Tests: toroidal surfaces, cylindrical surfaces
   - Export at: Fine mesh

6. **bracket_simple.stl** - L-bracket with holes
   - Tests: multiple planes at angles, cylindrical holes
   - Export at: Normal mesh

### Varying Triangle Sizes
7. **plate_with_hole_100x50.stl** - 100x50mm plate with 10mm center hole
   - Tests: large flat faces vs small curved facesc
   - Important: Do NOT use adaptive mesh - want to see large triangles on plate
   - Export at: Coarse mesh

8. **stepped_block.stl** - Block with steps of varying sizes
   - Tests: many parallel planes at different heights
   - Export at: Coarse mesh (want clean planar triangles)

### Complex Surfaces
9. **dome_hemisphere_20mm.stl** - 20mm diameter hemisphere
   - Tests: large curved surface recognition
   - Export at: Fine mesh

10. **chamfered_cube_10mm_c1.stl** - 10mm cube with 1mm edge chamfers
    - Tests: angled plane detection, small faces
    - Export at: Normal mesh

## Export Settings

### Fusion 360
- Format: STL (ASCII or Binary)
- Refinement: Choose based on test needs:
  - Coarse: For testing large face handling
  - Medium: General testing
  - Fine: For testing curved surface detail
- Units: Millimeters (default)

### Onshape
- Format: STL
- Units: Millimeters
- Resolution:
  - Coarse: For flat surfaces
  - Medium: General testing  
  - Fine: For curved surfaces
- Triangle count: Note the count for expected behavior

## Naming Convention

```
<shape>_<dimension_mm>[_<feature>].stl
```

Examples:
- `sphere_25mm.stl`
- `cylinder_10x30mm.stl` (diameter x height)
- `cube_10mm_r2_fillet.stl` (10mm cube, 2mm fillet radius)

## Adding New Test Files

When adding new test STL files:
1. Place in appropriate subdirectory (`manual/` or `generated/`)
2. Update this README with:
   - Filename and dimensions
   - What surface types it tests
   - Export settings used
3. Consider adding a corresponding unit test if testing new functionality
