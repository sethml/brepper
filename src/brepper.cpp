#include "brepper.hpp"
#include "common/logging.hpp"
#include "mesh/stl_reader.hpp"
#include "mesh/mesh_sampling.hpp"
#include "segmentation/ransac_segmenter.hpp"
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/conversions.h>
#include <iostream>
#include <iomanip>
#include <limits>

namespace brepper {

BrepperPipeline::BrepperPipeline(const Config& config) : config_(config) {
    // Set up logging based on config
    Logger& logger = Logger::instance();
    logger.set_quiet(config.quiet);
    
    if (config.debug) {
        logger.set_level(LogLevel::DEBUG);
    } else if (config.verbose) {
        logger.set_level(LogLevel::INFO);
    } else {
        logger.set_level(LogLevel::WARNING);
    }
}

bool BrepperPipeline::process() {
    LOG_INFO("Starting brepper pipeline");
    LOG_INFO("Input: ", config_.input_file);
    if (config_.stop_after_stage == PipelineStage::Export) {
        LOG_INFO("Output: ", config_.output_file);
    } else {
        LOG_INFO("Stopping after stage ", static_cast<int>(config_.stop_after_stage));
    }
    
    // Execute pipeline stages up to the requested stage
    if (!stage1_load_mesh()) {
        LOG_ERROR("Stage 1 (load mesh) failed");
        return false;
    }
    
    // Print dimensions if requested
    if (config_.print_dimensions) {
        print_mesh_dimensions();
    }
    
    if (config_.stop_after_stage < PipelineStage::Segment) {
        LOG_INFO("Pipeline completed (stopped after stage 1)");
        return true;
    }
    
    if (!stage2_segment_surfaces()) {
        LOG_ERROR("Stage 2 (segment surfaces) failed");
        return false;
    }
    
    if (config_.stop_after_stage < PipelineStage::Assign) {
        LOG_INFO("Pipeline completed (stopped after stage 2)");
        return true;
    }
    
    if (!stage3_assign_triangles()) {
        LOG_ERROR("Stage 3 (assign triangles) failed");
        return false;
    }
    
    if (config_.stop_after_stage < PipelineStage::Boundary) {
        LOG_INFO("Pipeline completed (stopped after stage 3)");
        return true;
    }
    
    if (!stage4_detect_boundaries()) {
        LOG_ERROR("Stage 4 (detect boundaries) failed");
        return false;
    }
    
    if (config_.stop_after_stage < PipelineStage::BRep) {
        LOG_INFO("Pipeline completed (stopped after stage 4)");
        return true;
    }
    
    if (!stage5_build_brep()) {
        LOG_ERROR("Stage 5 (build B-Rep) failed");
        return false;
    }
    
    if (config_.stop_after_stage < PipelineStage::Export) {
        LOG_INFO("Pipeline completed (stopped after stage 5)");
        return true;
    }
    
    if (!stage6_export_step()) {
        LOG_ERROR("Stage 6 (export STEP) failed");
        return false;
    }
    
    LOG_INFO("Pipeline completed successfully");
    return true;
}

bool BrepperPipeline::stage1_load_mesh() {
    LOG_INFO("Stage 1: Loading and preprocessing mesh");
    
    // Step 1.1: Load STL file
    STLReader reader(config_);
    if (!reader.load(config_.input_file, results_.input_mesh)) {
        return false;
    }
    
    // Step 1.2: Sample mesh to point cloud with normals
    MeshSampler sampler(config_);
    if (!sampler.sample(results_.input_mesh, results_.sampled_cloud)) {
        return false;
    }
    
    // Step 1.3: Save debug output if requested
    if (!config_.save_point_cloud.empty()) {
        LOG_INFO("Saving point cloud to: ", config_.save_point_cloud);
        
        // Determine format from extension
        std::string ext = config_.save_point_cloud.substr(
            config_.save_point_cloud.find_last_of('.') + 1);
        
        if (ext == "pcd") {
            pcl::io::savePCDFileBinary(config_.save_point_cloud, *results_.sampled_cloud);
        } else if (ext == "ply") {
            pcl::io::savePLYFileBinary(config_.save_point_cloud, *results_.sampled_cloud);
        } else {
            LOG_WARN("Unknown point cloud format, saving as PCD");
            pcl::io::savePCDFileBinary(config_.save_point_cloud + ".pcd", *results_.sampled_cloud);
        }
    }
    
    LOG_INFO("Stage 1 complete: ", results_.sampled_cloud->size(), " points generated");
    return true;
}

bool BrepperPipeline::stage2_segment_surfaces() {
    LOG_INFO("Stage 2: Segmenting surfaces with RANSAC");
    
    RANSACSegmenter segmenter(config_);
    if (!segmenter.segment(results_.sampled_cloud, results_.fitted_surfaces)) {
        LOG_ERROR("RANSAC segmentation failed");
        return false;
    }
    
    LOG_INFO("Stage 2 complete: found ", results_.fitted_surfaces.size(), " surfaces");
    return true;
}

bool BrepperPipeline::stage3_assign_triangles() {
    LOG_INFO("Stage 3: Assigning triangles to surfaces");
    // TODO: Implement triangle assignment
    LOG_WARN("Stage 3: Not implemented yet");
    return true;
}

bool BrepperPipeline::stage4_detect_boundaries() {
    LOG_INFO("Stage 4: Detecting boundaries and fitting curves");
    // TODO: Implement boundary detection
    LOG_WARN("Stage 4: Not implemented yet");
    return true;
}

bool BrepperPipeline::stage5_build_brep() {
    LOG_INFO("Stage 5: Building B-Rep from surfaces");
    // TODO: Implement OCCT B-Rep construction
    LOG_WARN("Stage 5: Not implemented yet");
    return true;
}

bool BrepperPipeline::stage6_export_step() {
    LOG_INFO("Stage 6: Exporting STEP file");
    // TODO: Implement STEP export
    LOG_WARN("Stage 6: Not implemented yet");
    return true;
}

void BrepperPipeline::print_mesh_dimensions() {
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromPCLPointCloud2(results_.input_mesh.cloud, *vertices);
    
    float min_x = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::lowest();
    float min_y = std::numeric_limits<float>::max();
    float max_y = std::numeric_limits<float>::lowest();
    float min_z = std::numeric_limits<float>::max();
    float max_z = std::numeric_limits<float>::lowest();
    
    for (const auto& pt : *vertices) {
        min_x = std::min(min_x, pt.x); max_x = std::max(max_x, pt.x);
        min_y = std::min(min_y, pt.y); max_y = std::max(max_y, pt.y);
        min_z = std::min(min_z, pt.z); max_z = std::max(max_z, pt.z);
    }
    
    // Note: coordinates are already converted to mm during load
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Mesh dimensions (mm):\n";
    std::cout << "  X: " << (max_x - min_x) << " mm  (range: " << min_x << " to " << max_x << ")\n";
    std::cout << "  Y: " << (max_y - min_y) << " mm  (range: " << min_y << " to " << max_y << ")\n";
    std::cout << "  Z: " << (max_z - min_z) << " mm  (range: " << min_z << " to " << max_z << ")\n";
}

} // namespace brepper