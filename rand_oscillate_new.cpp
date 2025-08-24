// rand_oscillate_new.cpp
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#if defined(__APPLE__)
#  include <mach/mach.h>
#  include <malloc/malloc.h>
#  include <unistd.h>
#elif defined(__linux__)
#  include <unistd.h>
#  include <sys/types.h>
#  include <sys/sysinfo.h>
#endif

// 返回更接近真实占用的指标（macOS: phys_footprint；Linux: RSS）
static inline std::size_t get_rss_like_bytes() {
#if defined(__APPLE__)
    task_vm_info_data_t info{};
    mach_msg_type_number_t cnt = TASK_VM_INFO_COUNT;
    if (task_info(mach_task_self(), TASK_VM_INFO, (task_info_t)&info, &cnt) == KERN_SUCCESS) {
        return (std::size_t)info.phys_footprint; // 推荐用这个
    }
    return 0;
#elif defined(__linux__)
    long rss_pages = 0;
    FILE *f = std::fopen("/proc/self/statm", "r");
    if (f) {
        long total_pages = 0;
        if (std::fscanf(f, "%ld %ld", &total_pages, &rss_pages) != 2) rss_pages = 0;
        std::fclose(f);
    }
    return (std::size_t)rss_pages * (std::size_t)sysconf(_SC_PAGESIZE);
#else
    return 0;
#endif
}

static inline const char* human(std::size_t bytes, double& v) {
    static const char* unit[] = { "B","KB","MB","GB","TB" };
    v = (double)bytes;
    int i = 0;
    while (v >= 1024.0 && i < 4) { v /= 1024.0; ++i; }
    return unit[i];
}

struct Block {
    std::unique_ptr<std::uint8_t[]> ptr;
    std::size_t bytes{};
};

int main(int argc, char** argv) {
    // 参数: iterations minMiB maxMiB sleep_ms relief_flag
    int iterations  = (argc > 1) ? std::atoi(argv[1]) : 40;
    std::size_t minMiB = (argc > 2) ? std::strtoull(argv[2], nullptr, 10) : 4;
    std::size_t maxMiB = (argc > 3) ? std::strtoull(argv[3], nullptr, 10) : 1024;
    int sleep_ms   = (argc > 4) ? std::atoi(argv[4]) : 800;
    bool use_relief = (argc > 5) ? (std::atoi(argv[5]) != 0) : true;

#if defined(__APPLE__) || defined(__linux__)
    std::size_t page = (std::size_t)getpagesize();
#else
    std::size_t page = 4096;
#endif

    std::mt19937_64 rng{ std::random_device{}() };
    // 用对数均匀分布（更常抽到小块，偶尔抽到超大块，更贴近真实 workload）
    std::uniform_real_distribution<double> U( std::log((double)minMiB), std::log((double)maxMiB) );

    std::vector<Block> live; live.reserve(128);

    std::puts("Random allocate & occasionally free the current largest block.");
    std::printf("iterations=%d, range=[%zu,%zu] MiB, sleep=%dms, pressure_relief=%s\n",
                iterations, minMiB, maxMiB, sleep_ms, use_relief ? "ON" : "OFF");
    std::puts("Tips: try `MallocNanoZone=0 MallocSpaceEfficient=1` to see clearer drops.");

    for (int it = 1; it <= iterations; ++it) {
        // --- 随机大小（log-uniform 到 MiB） ---
        std::size_t szMiB = (std::size_t)std::llround(std::exp(U(rng)));
        szMiB = std::max<std::size_t>(szMiB, 1);
        std::size_t bytes = szMiB * 1024ull * 1024ull;

        // --- 分配并触页 ---
        Block b;
        b.ptr.reset(new std::uint8_t[bytes]);
        b.bytes = bytes;
        for (std::size_t off = 0; off < bytes; off += page) b.ptr[off] = 1;

        double v; const char* u;
        u = human(get_rss_like_bytes(), v);
        std::printf("[iter %2d] alloc %6zu MiB -> RSS %.2f %s (live=%zu)\n",
                    it, szMiB, v, u, live.size()+1);
        std::fflush(stdout);

        live.emplace_back(std::move(b));

        // --- 每隔 3 轮，释放当前最大的块，观察回落 ---
        if (it % 3 == 0 && !live.empty()) {
            // 找最大的
            auto itMax = std::max_element(live.begin(), live.end(),
                [](const Block& a, const Block& b){ return a.bytes < b.bytes; });
            std::size_t freeMiB = itMax->bytes / (1024ull*1024ull);

            u = human(get_rss_like_bytes(), v);
            std::printf("          before free largest(%zu MiB): RSS %.2f %s\n", freeMiB, v, u);

            // delete[]
            itMax->ptr.reset();
            // 从 live 中移除条目（保持容器干净）
            live.erase(itMax);

            u = human(get_rss_like_bytes(), v);
            std::printf("          after  delete[] largest:      RSS %.2f %s\n", v, u);

#if defined(__APPLE__)
            if (use_relief) {
                malloc_zone_pressure_relief(nullptr, 0);
                u = human(get_rss_like_bytes(), v);
                std::printf("          after  pressure_relief():     RSS %.2f %s\n", v, u);
            }
#endif
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }

    // 结尾把剩余全部释放并观测
    double v; const char* u;
    u = human(get_rss_like_bytes(), v);
    std::printf("== before final clear: RSS %.2f %s (live=%zu)\n", v, u, live.size());
    live.clear();
#if defined(__APPLE__)
    if (use_relief) {
        malloc_zone_pressure_relief(nullptr, 0);
    }
#endif
    u = human(get_rss_like_bytes(), v);
    std::printf("== after  final clear: RSS %.2f %s (live=%zu)\n", v, u, live.size());
}