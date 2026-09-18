#pragma once

// Minimal replacement for D's std.experimental.logger.  Log lines are written
// to a file (created along with its parent directories) in the same style as
// the original: "2020-05-25T16:02:22.223 [info] message...".
//
// fatal() is [[noreturn]]: it logs the message and aborts the process.  This
// mirrors D's std.experimental.logger.fatal, which also terminates.  Note that
// RAII cleanup (e.g. nc_close guards) intentionally does NOT run on fatal
// paths -- the process is dying anyway.

#include <array>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace niceops {

namespace detail {

// Streaming helpers so log calls can take containers (D's logger prints arrays
// inline).  The vector/array overloads recurse into the elements.
template <typename T>
void print(std::ostream& os, const T& value)
{
    os << value;
}

template <typename T>
void print(std::ostream& os, const std::vector<T>& values)
{
    os << "[";
    for (std::size_t i = 0; i < values.size(); i++)
    {
        if (i)
            os << ", ";
        print(os, values[i]);
    }
    os << "]";
}

template <typename T, std::size_t N>
void print(std::ostream& os, const std::array<T, N>& values)
{
    os << "[";
    for (std::size_t i = 0; i < N; i++)
    {
        if (i)
            os << ", ";
        print(os, values[i]);
    }
    os << "]";
}

} // namespace detail

enum class LogLevel
{
    trace = 0,
    info,
    warn,
    error,
    fatal
};

class Logger
{
public:
    static Logger& instance()
    {
        static Logger logger;
        return logger;
    }

    void set_file(const std::string& path)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (file_.is_open())
            file_.close();
        auto parent = std::filesystem::path(path).parent_path();
        if (!parent.empty())
            std::filesystem::create_directories(parent);
        file_.open(path, std::ios::app);
    }

    void set_level(LogLevel level)
    {
        min_level_ = level;
    }

    template <typename... Args>
    void write(LogLevel level, Args&&... args)
    {
        if (level < min_level_)
            return;
        std::ostringstream oss;
        join(oss, std::forward<Args>(args)...);
        std::lock_guard<std::mutex> lock(mutex_);
        auto line = timestamp() + " [" + level_name(level) + "] " + oss.str() + "\n";
        if (file_.is_open())
        {
            file_ << line;
            file_.flush();
        }
        else
        {
            std::cerr << line;
        }
    }

private:
    Logger() = default;

    static const char* level_name(LogLevel level)
    {
        switch (level)
        {
        case LogLevel::trace:
            return "trace";
        case LogLevel::info:
            return "info";
        case LogLevel::warn:
            return "warn";
        case LogLevel::error:
            return "error";
        case LogLevel::fatal:
            return "fatal";
        }
        return "?";
    }

    static std::string timestamp()
    {
        auto now = std::chrono::system_clock::now();
        auto t = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
        std::tm tm {};
#if defined(_WIN32)
        gmtime_s(&tm, &t);
#else
        gmtime_r(&t, &tm);
#endif
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
        char msbuf[8];
        std::snprintf(msbuf, sizeof(msbuf), ".%03d", static_cast<int>(ms.count()));
        return std::string(buf) + msbuf;
    }

    void join(std::ostringstream&) {}
    template <typename T, typename... Rest>
    void join(std::ostringstream& oss, T&& value, Rest&&... rest)
    {
        detail::print(oss, value);
        join(oss, std::forward<Rest>(rest)...);
    }

    std::ofstream file_;
    std::mutex mutex_;
    LogLevel min_level_ = LogLevel::info;
};

template <typename... Args>
void trace(Args&&... args)
{
    Logger::instance().write(LogLevel::trace, std::forward<Args>(args)...);
}

// NOTE: log() and info() are the same level.  Both names are in use across the
// codebase; they are aliases, not a severity distinction.
template <typename... Args>
void log(Args&&... args)
{
    Logger::instance().write(LogLevel::info, std::forward<Args>(args)...);
}

template <typename... Args>
void info(Args&&... args)
{
    Logger::instance().write(LogLevel::info, std::forward<Args>(args)...);
}

template <typename... Args>
void warn(Args&&... args)
{
    Logger::instance().write(LogLevel::warn, std::forward<Args>(args)...);
}

template <typename... Args>
void error(Args&&... args)
{
    Logger::instance().write(LogLevel::error, std::forward<Args>(args)...);
}

template <typename... Args>
[[noreturn]] void fatal(Args&&... args)
{
    Logger::instance().write(LogLevel::fatal, std::forward<Args>(args)...);
    std::abort();
}

} // namespace niceops
