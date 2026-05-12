#ifndef __MYTIMER_H_
#define __MYTIMER_H_

#include <iostream>
#include <chrono>
#include <string>
#include <iomanip>
#include <sstream>
#include <ctime>

class Timer
{
private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    TimePoint start_time_;
    bool is_running_;

    // Convert seconds to human-readable time string
    std::string formatTime(double seconds) const
    {
        int days = static_cast<int>(seconds) / 86400;
        int hours = (static_cast<int>(seconds) % 86400) / 3600;
        int minutes = (static_cast<int>(seconds) % 3600) / 60;
        int secs = static_cast<int>(seconds) % 60;

        std::stringstream ss;
        if (days > 0)
        {
            ss << days << "d";
            if (hours > 0)
            {
                ss << " " << hours << "h";
            }
        }
        else if (hours > 0)
        {
            ss << hours << "h";
            if (minutes > 0)
            {
                ss << " " << minutes << "m";
            }
        }
        else if (minutes > 0)
        {
            ss << minutes << "m";
            if (secs > 0)
            {
                ss << " " << secs << "s";
            }
        }
        else
        {
            ss << secs << "s";
        }
        return ss.str();
    }

    // Get formatted completion time
    std::string getCompletionTime(double eta_seconds) const
    {
        auto now = std::chrono::system_clock::now();
        auto eta_duration = std::chrono::seconds(static_cast<int>(eta_seconds));
        auto completion_time = now + std::chrono::duration_cast<std::chrono::system_clock::duration>(eta_duration);

        std::time_t completion_time_t = std::chrono::system_clock::to_time_t(completion_time);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&completion_time_t), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }

public:
    // Constructor
    Timer() : is_running_(false) {}

    // Start the timer
    void start()
    {
        start_time_ = Clock::now();
        is_running_ = true;
    }

    // Stop the timer and return elapsed time in seconds
    double stop()
    {
        if (is_running_)
        {
            Duration elapsed = Clock::now() - start_time_;
            is_running_ = false;
            return elapsed.count();
        }
        return 0.0;
    }

    // Get elapsed time without stopping the timer
    double elapsed() const
    {
        if (is_running_)
        {
            Duration elapsed = Clock::now() - start_time_;
            return elapsed.count();
        }
        return 0.0;
    }

    // Print elapsed time with message
    void printElapsed() const
    {
        double elapsed_time = elapsed();
        double eta = (elapsed_time > 0) ? (elapsed_time / elapsed_time) * (1 - elapsed_time) : 0;
        std::string elapsed_str = formatTime(elapsed_time);
        std::string eta_str = formatTime(eta);

        // Get completion time
        std::string completion_time = getCompletionTime(eta);

        std::cout << std::fixed << std::setprecision(3) << (elapsed_time / elapsed_time) * 100.0 << "% "
                  << "Elapsed: " << elapsed_str << " "
                  << "ETA: " << eta_str << " "
                  << "Complete at: " << completion_time
                  << std::endl;
    }

    // Update progress by current progress fraction (0 to 1)
    void updateProgress(double progress)
    {
        if (is_running_)
        {
            double elapsed_time = elapsed();
            double eta = (progress > 0) ? (elapsed_time / progress) * (1 - progress) : 0;
            std::string elapsed_str = formatTime(elapsed_time);
            std::string eta_str = formatTime(eta);

            // Get completion time
            std::string completion_time = getCompletionTime(eta);

            std::cout << "\r" // 返回行首，准备覆盖之前的输出
                      << "Progress: \033[1;32m" << std::fixed << std::setprecision(2) << progress * 100.0 << "%\033[0m "
                      << "| Elapsed: \033[1;34m" << elapsed_str << "\033[0m "
                      << "| ETA: \033[1;33m" << eta_str << "\033[0m "
                      << "| Complete at: \033[1;35m" << completion_time << "\033[0m"
                      << std::flush;
        }
    }
};

#endif
