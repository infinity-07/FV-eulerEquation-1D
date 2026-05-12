#ifndef _MYTIME_H_
#define _MYTIME_H_

#include <iostream>
#include <ctime>
#include <sstream>

class CTime
{

private:
    std::clock_t m_startTime;   // 计时器的起始时间
    std::clock_t m_currentTime; // 当前时间
    std::clock_t m_finishTime;  // 计时器的结束时间
    bool m_isStarted = false;   // 计时器的状态标记，记录计时器是否已经开始

public:
    double m_totalSeconds;       // 累计的总时间（秒）
    std::string m_startTimeStr;  // 计时器开始时间的字符串表示
    std::string m_finishTimeStr; // 计时器结束时间的字符串表示

public:
    // 开始计时
    void start()
    {
        if (m_isStarted)
        {
            std::cout << "Timer is already started." << std::endl; // 如果计时器已经开始，则输出提示信息
            return;
        }

        m_startTime = std::clock(); // 获取当前时间作为计时器的起始时间
        m_isStarted = true;         // 将计时器状态标记为已开始

        // 获取当前时间
        std::time_t now = std::time(0);

        // 将时间转换为字符串
        m_startTimeStr = getTimeString(now);
    }

    // 获取当前已累计的时间（秒）
    double getCurrent()
    {
        if (!m_isStarted)
        {
            std::cout << "Timer is not started." << std::endl; // 如果计时器未开始，则输出提示信息
            return 0.0;
        }

        m_currentTime = std::clock(); // 获取当前时间

        return double(m_currentTime - m_startTime) / CLOCKS_PER_SEC; // 返回时间差（秒）
    }

    // 终止计时器，输出程序运行消耗的时间
    void end()
    {
        if (!m_isStarted)
        {
            std::cout << "Timer is not started." << std::endl; // 如果计时器未开始，则输出提示信息
            return;
        }

        m_finishTime = std::clock(); // 获取当前时间作为计时器的结束时间
        m_isStarted = false;         // 将计时器状态标记为未开始

        m_totalSeconds = double(m_finishTime - m_startTime) / CLOCKS_PER_SEC;

        // 获取当前时间
        std::time_t now = std::time(0);

        // 将时间转换为字符串
        m_finishTimeStr = getTimeString(now);
    }

    // // 输出计时器的运行时间
    // void output()
    // {
    //     std::cout << "========================================" << std::endl;
    //     std::cout << "Running time information:" << std::endl;
    //     std::cout << "Total time is: " << m_totalSeconds << " s" << std::endl;             // 输出总时间（秒）
    //     std::cout << "Total time is: " << m_totalSeconds / 60.0 << " min" << std::endl;    // 输出总时间（分钟）
    //     std::cout << "Total time is: " << m_totalSeconds / 3600.0 << " hour" << std::endl; // 输出总时间（小时）
    //     std::cout << "========================================" << std::endl;
    //     std::cout << std::endl;
    // }

    // 将 std::time_t 类型的时间转换为字符串表示
    std::string getTimeString(const std::time_t &time)
    {
        std::stringstream ss;
        tm *timeInfo = std::localtime(&time);

        ss << timeInfo->tm_year + 1900 << "年"
           << timeInfo->tm_mon + 1 << "月"
           << timeInfo->tm_mday << "日"
           << timeInfo->tm_hour << ":"
           << timeInfo->tm_min << ":"
           << timeInfo->tm_sec;

        return ss.str();
    }

    // 将秒数转换为带单位的字符串表示
    std::string getTotalTimeString(double secondTime)
    {
        std::string result;
        std::string unit;

        if (secondTime < 60)
        {
            result = std::to_string(secondTime);
            unit = "s";
        }
        else if (secondTime < 3600)
        {
            double minutes = secondTime / 60;
            result = std::to_string(minutes);
            unit = "分钟";
        }
        else if (secondTime < 3600 * 24)
        {
            double hours = secondTime / 3600;
            result = std::to_string(hours);
            unit = "小时";
        }
        else
        {
            double days = secondTime / (3600 * 24);
            result = std::to_string(days);
            unit = "天";
        }

        // 保留两位小数
        size_t decimalPos = result.find(".");
        if (decimalPos != std::string::npos && decimalPos + 3 < result.length())
        {
            result = result.substr(0, decimalPos + 3);
        }

        return result + unit;
    }
};

#endif
