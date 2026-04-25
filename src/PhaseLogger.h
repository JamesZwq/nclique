// PhaseLogger.h
//
// Lightweight per-phase time + RSS + component-bytes logger for the breakdown
// experiment.  Header-only; one TU's worth of static state lives in the
// inline singleton.  Activated by setting the env var PIVOTER_BREAKDOWN_LOG
// to a file path; if unset, the logger is essentially free (still snapshots
// time/RSS at each mark, but never writes anything).
//
// Usage:
//   #include "PhaseLogger.h"
//   ...
//   daf::phaseStart();
//   ... work ...
//   daf::phaseMark("loadAndSort");
//   ... more work ...
//   daf::phaseMark("buildSDCT", treeBytes);  // optional component byte tag
//   ... peel ...
//   daf::phaseMark("peel");
//   daf::phaseDump();    // appends one TSV row per recorded phase to the log
//
// The TSV header (written once if file is empty) is:
//   meta\tphase\tduration_ms\trss_kb\tdelta_rss_kb\tcomponent_bytes
//
// The `meta` column is taken verbatim from PIVOTER_BREAKDOWN_META, so the
// runner script encodes (graph, r, s, algo, run_id, ...) into a single
// comma-separated tag.
//

#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#if defined(__APPLE__)
#  include <mach/mach.h>
#endif

namespace daf {

inline long phaseLogger_currentRSS_kB() {
#if defined(__linux__)
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            long kb = 0;
            std::istringstream iss(line.substr(6));
            iss >> kb;
            return kb;
        }
    }
    return 0;
#elif defined(__APPLE__)
    task_basic_info t_info;
    mach_msg_type_number_t cnt = TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_BASIC_INFO,
                  (task_info_t)&t_info, &cnt) == KERN_SUCCESS) {
        return (long)(t_info.resident_size / 1024);
    }
    return 0;
#else
    return 0;
#endif
}

struct PhaseRecord {
    std::string name;
    double duration_ms;
    long rss_kb;
    long delta_rss_kb;
    long component_bytes;
};

class PhaseLogger {
public:
    static PhaseLogger& instance() {
        static PhaseLogger inst;
        return inst;
    }

    void start() {
        t_prev_  = std::chrono::steady_clock::now();
        rss_prev_ = phaseLogger_currentRSS_kB();
        records_.clear();
        active_ = true;
    }

    void mark(const std::string& name, long component_bytes = 0) {
        if (!active_) return;
        auto now = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(now - t_prev_).count();
        const long rss = phaseLogger_currentRSS_kB();
        records_.push_back({name, ms, rss, rss - rss_prev_, component_bytes});
        t_prev_  = now;
        rss_prev_ = rss;
    }

    void dump() {
        if (!active_) return;
        const char* path = std::getenv("PIVOTER_BREAKDOWN_LOG");
        if (!path) return;
        const char* meta = std::getenv("PIVOTER_BREAKDOWN_META");
        // Open in append mode, write header if file is empty.
        bool needs_header = false;
        {
            std::ifstream probe(path, std::ios::ate);
            needs_header = !probe.good() || probe.tellg() == 0;
        }
        std::ofstream f(path, std::ios::app);
        if (!f) {
            std::fprintf(stderr, "[PhaseLogger] cannot open %s\n", path);
            return;
        }
        if (needs_header) {
            f << "meta\tphase\tduration_ms\trss_kb\tdelta_rss_kb\tcomponent_bytes\n";
        }
        for (const auto& r : records_) {
            f << (meta ? meta : "") << '\t'
              << r.name           << '\t'
              << r.duration_ms    << '\t'
              << r.rss_kb         << '\t'
              << r.delta_rss_kb   << '\t'
              << r.component_bytes << '\n';
        }
        records_.clear();
    }

private:
    PhaseLogger() = default;
    bool active_ = false;
    std::chrono::steady_clock::time_point t_prev_;
    long rss_prev_ = 0;
    std::vector<PhaseRecord> records_;
};

inline void phaseStart()                                   { PhaseLogger::instance().start(); }
inline void phaseMark(const std::string& name, long bytes = 0) {
    PhaseLogger::instance().mark(name, bytes);
}
inline void phaseDump()                                    { PhaseLogger::instance().dump(); }

} // namespace daf
