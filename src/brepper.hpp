#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

namespace brepper {

class BrepperPipeline {
public:
    explicit BrepperPipeline(const Config& config);
    
    // Main processing pipeline
    bool process();
    
    // Access results
    const ProcessingResults& results() const { return results_; }

private:
    // Pipeline stages
    bool stage1_load_mesh();
    bool stage2_segment_surfaces(); 
    bool stage3_assign_triangles();
    bool stage4_detect_boundaries();
    bool stage5_build_brep();
    bool stage6_export_step();
    
    // Utilities
    void print_mesh_dimensions();
    
    Config config_;
    ProcessingResults results_;
};

} // namespace brepper