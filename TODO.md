# Ports from C++ Implementation

## Full pipeline tests for cone and torus models
All test models including cone and torus models and filleted_cylinder, filleted_hole_block, filleted_pipe pass --compare through stage 4.1. Unblock existing `pipe_elbow_10_fine` (onshape) test.

# Unified Surface-of-Revolution Framework (refactor)

## Extract shared SOR axis estimation
Extract shared SOR axis estimation (normal covariance smallest eigenvector) into a common utility function.

## Implement (h,r) coordinate reduction
For a set of vertices and an axis, compute axial and radial coordinates.

## Implement unified profile classifier
Given (h,r) points, fit horizontal line (cylinder), sloped line (cone), circle centered on axis (sphere), and offset circle (torus). Return all fits with RMS residuals.

## Refactor cylinder fitting to use (h,r) reduction
Use (h,r) reduction + profile classification as an alternative initialization.

## Refactor sphere fitting to use (h,r) reduction
Use (h,r) reduction when axis is estimated from normal covariance.

## Integrate unified SOR into BFS
After initial seed, try all profile types and pick best before growing.

## Regression tests for SOR refactor
Ensure all existing test models still pass after refactor.

# Gaussian Map Pre-Classification

## Implement normal covariance eigenvalue classification
Compute λ₁, λ₂, λ₃ eigenvalue pattern for face neighborhoods. Classify as planar (λ₁ ≫ λ₂ ≈ λ₃), cylindrical/conical (λ₁ ≈ λ₂ ≫ λ₃), spherical (λ₁ ≈ λ₂ ≈ λ₃), or toroidal.

## Use classification to prioritize seeding strategies
Prioritize seeding strategies per region, reducing wasted fitting attempts.

## Add quadric fitting as fast initialization
10-coefficient SVD as a fast initialization alternative. Extract surface type from eigenvalue pattern of the 3×3 quadric matrix.

# NURBS Fallback

## Create test models with freeform surfaces
Create test models with freeform surfaces that cannot be fit by analytic types (e.g., sculpted surface, spline loft).

## Implement NURBS fitting
Implement NURBS fitting for remaining unassigned face groups using OCCT's `GeomAPI_PointsToBSplineSurface`.

## Integrate NURBS into stage 3 reconstruction
Integrate NURBS surfaces into stage 3 reconstruction (edge curves, face creation).

## Extend surface selection for NURBS
Extend stage 2.6 surface selection to handle NURBS alongside analytic surfaces.

# Future Enhancements

## Support for OBJ, PLY input formats

## IGES export option
