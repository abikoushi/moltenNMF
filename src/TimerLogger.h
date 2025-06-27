#include <chrono>
#include <iostream>
#include <fstream>
#include <mutex>

class TimerLogger {
public:
  TimerLogger(const std::string& label, const std::string& filename = "timelog.csv")
    : label_(label), filename_(filename) {
    start_ = std::chrono::high_resolution_clock::now();
  }
  
  ~TimerLogger() {
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count() / 1000.0;
    log_to_csv(label_, elapsed);
  }
  
private:
  std::string label_;
  std::string filename_;
  std::chrono::high_resolution_clock::time_point start_;
  
  static void log_to_csv(const std::string& label, double elapsed) {
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    
    std::ofstream ofs("timelog.csv", std::ios::app);
    if (ofs.is_open()) {
      ofs << label << "," << elapsed << "\n";
    }
  }
};
