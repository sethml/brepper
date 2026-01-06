#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace brepper {

enum class LogLevel {
    DEBUG = 0,
    INFO = 1, 
    WARNING = 2,
    ERROR = 3
};

class Logger {
public:
    static Logger& instance() {
        static Logger logger;
        return logger;
    }
    
    void set_level(LogLevel level) { 
        current_level_ = level; 
    }
    
    void set_quiet(bool quiet) { 
        quiet_ = quiet; 
    }
    
    template<typename... Args>
    void debug(Args&&... args) {
        log(LogLevel::DEBUG, "[DEBUG] ", std::forward<Args>(args)...);
    }
    
    template<typename... Args>
    void info(Args&&... args) {
        log(LogLevel::INFO, "[INFO]  ", std::forward<Args>(args)...);
    }
    
    template<typename... Args>
    void warn(Args&&... args) {
        log(LogLevel::WARNING, "[WARN]  ", std::forward<Args>(args)...);
    }
    
    template<typename... Args>
    void error(Args&&... args) {
        log(LogLevel::ERROR, "[ERROR] ", std::forward<Args>(args)...);
    }

private:
    LogLevel current_level_ = LogLevel::INFO;
    bool quiet_ = false;
    
    template<typename... Args>
    void log(LogLevel level, const std::string& prefix, Args&&... args) {
        if (quiet_ && level < LogLevel::ERROR) return;
        if (level < current_level_) return;
        
        std::ostringstream oss;
        oss << prefix;
        (oss << ... << args);
        oss << std::endl;
        
        if (level == LogLevel::ERROR) {
            std::cerr << oss.str();
        } else {
            std::cout << oss.str();
        }
    }
};

// Convenience macros
#define LOG_DEBUG(...) brepper::Logger::instance().debug(__VA_ARGS__)
#define LOG_INFO(...) brepper::Logger::instance().info(__VA_ARGS__)  
#define LOG_WARN(...) brepper::Logger::instance().warn(__VA_ARGS__)
#define LOG_ERROR(...) brepper::Logger::instance().error(__VA_ARGS__)

} // namespace brepper