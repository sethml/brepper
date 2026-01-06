#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

namespace brepper {

class RANSACSegmenter {
public:
    explicit RANSACSegmenter(const Config& config);
    
    // Segment surfaces from point cloud using iterative RANSAC
    bool segment(PointCloudNormalPtr cloud, std::vector<FittedSurface>& surfaces);
    
private:
    const Config& config_;
    
    // Fit single surface type and return best result
    bool fit_surface_type(
        PointCloudNormalPtr cloud,
        SurfaceType type,
        FittedSurface& surface
    );
    
    bool fit_plane(PointCloudNormalPtr cloud, FittedSurface& surface);
    bool fit_cylinder(PointCloudNormalPtr cloud, FittedSurface& surface);
    bool fit_sphere(PointCloudNormalPtr cloud, FittedSurface& surface);
    bool fit_cone(PointCloudNormalPtr cloud, FittedSurface& surface);
};

} // namespace brepper