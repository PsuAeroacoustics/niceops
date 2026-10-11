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
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "paths.hpp"

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

// Build one message out of the argument pack, using the printers above.
inline void join(std::ostringstream&) {}

template <typename T, typename... Rest>
void join(std::ostringstream& oss, T&& value, Rest&&... rest)
{
    print(oss, value);
    join(oss, std::forward<Rest>(rest)...);
}

template <typename... Args>
std::string format_message(Args&&... args)
{
    std::ostringstream oss;
    join(oss, std::forward<Args>(args)...);
    return oss.str();
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

namespace detail {

// The active minimum level, mirrored out of the Logger so the free functions
// below can gate on it without touching the singleton.  Logger::instance() is
// a function-local static, so every call to it reloads and tests a guard
// variable -- affordable once, but trace() runs hundreds of millions of times
// inside the eval loop, where that guard and the argument formatting it
// protects were together costing several percent of the run.
inline LogLevel min_level = LogLevel::info;

} // namespace detail

// Whether a message at this level would be written at all.  Cheap enough to
// put in front of a hot-path log call, and it lets the compiler drop the
// argument evaluation when logging is off.
inline bool log_enabled(LogLevel level)
{
    return level >= detail::min_level;
}

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
        ensure_parent_directory(path);
        file_.open(path, std::ios::app);
    }

    // Back to stderr, as before any set_file.
    void close_file()
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (file_.is_open())
            file_.close();
    }

    void set_level(LogLevel level)
    {
        detail::min_level = level;
    }

    // Deliberately not inlined.  The formatting machinery below is large, and
    // when it was inlined into trace() the compiler stopped inlining trace()
    // itself -- leaving a real call, with its argument setup, at every hot-path
    // trace site.  Out of line, trace() shrinks to the level test, which does
    // fold into the caller.
    template <typename... Args>
#if defined(__GNUC__) || defined(__clang__)
    [[gnu::noinline]]
#endif
    void write(LogLevel level, Args&&... args)
    {
        // Warnings and errors are counted, and the first few kept, whatever -L
        // says, for the run manifest.
        if (level >= LogLevel::warn)
        {
            const auto message = detail::format_message(std::forward<Args>(args)...);
            note(level, message);
            write_line(level, message);
            return;
        }
        if (!log_enabled(level))
            return;
        write_line(level, detail::format_message(std::forward<Args>(args)...));
    }

    std::size_t warning_count() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return warnings_;
    }
    std::size_t error_count() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return errors_;
    }
    // The first warnings and errors, each prefixed with its level.
    std::vector<std::string> first_problems() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return first_problems_;
    }

    // The formatted half, once the message is a string.  input_error needs both
    // the log line and the same text for its exception.
    void write_line(LogLevel level, const std::string& message)
    {
        if (!log_enabled(level))
            return;
        std::lock_guard<std::mutex> lock(mutex_);
        auto line = timestamp() + " [" + level_name(level) + "] " + message + "\n";
        if (file_.is_open())
        {
            file_ << line;
            file_.flush();
            // The file is where a run's whole story goes, but nobody watches it
            // while the run is going: a skipped trajectory, a rejected option or
            // a fatal() abort used to leave the terminal empty and the exit
            // status the only sign anything was wrong.  Warnings and worse are
            // echoed to stderr as well; -L still decides what is written at all.
            if (level >= LogLevel::warn)
                std::cerr << line;
        }
        else
        {
            std::cerr << line;
        }
    }

private:
    Logger() = default;

    static constexpr std::size_t kept_problems = 20;

    void note(LogLevel level, const std::string& message)
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (level == LogLevel::warn)
            warnings_++;
        else
            errors_++;
        if (first_problems_.size() < kept_problems)
            first_problems_.push_back(std::string(level_name(level)) + ": " + message);
    }

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

    std::ofstream file_;
    mutable std::mutex mutex_;
    std::size_t warnings_ = 0;
    std::size_t errors_ = 0; // errors and fatals
    std::vector<std::string> first_problems_;
};

// Forced inline so the level test folds into the caller and the argument setup
// sinks into the branch that is not taken.  eval() and level_at_observer() make
// roughly seven trace calls each and run tens of millions of times, so leaving
// this as an ordinary call cost a sixth of the run doing nothing.
template <typename... Args>
#if defined(__GNUC__) || defined(__clang__)
[[gnu::always_inline]]
#endif
inline void trace(Args&&... args)
{
    if (!log_enabled(LogLevel::trace))
        return;
    Logger::instance().write(LogLevel::trace, std::forward<Args>(args)...);
}

template <typename... Args>
void info(Args&&... args)
{
    Logger::instance().write(LogLevel::info, std::forward<Args>(args)...);
}

// log() is info() under another name.  Both spellings are in use across the
// ported code and neither is worth churning, but they are one function rather
// than two identical ones, so a change to the level cannot reach only half the
// call sites.
template <typename... Args>
void log(Args&&... args)
{
    info(std::forward<Args>(args)...);
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

// What input_error throws: its message is already logged (and counted), so a
// caller reporting it again ("Skipping ...", "Stopping: ...") uses
// report_again rather than counting the one failure twice.
struct InputError : std::runtime_error
{
    using std::runtime_error::runtime_error;
};

// For a failure that condemns one input rather than the whole run: logs at
// error level and throws, so a caller processing several trajectories can note
// which one failed and carry on with the rest.  fatal() below still aborts, for
// the cases where continuing makes no sense.
template <typename... Args>
[[noreturn]] void input_error(Args&&... args)
{
    // Formatted once and used twice: the log line and the exception text cannot
    // drift apart.
    const auto message = detail::format_message(std::forward<Args>(args)...);
    Logger::instance().write(LogLevel::error, message);
    throw InputError(message);
}

// Logs a failure that is being reported a second time, with its context, at
// error level: counted only if it was not an InputError, which counted itself.
inline void report_again(const std::exception& problem, const std::string& context)
{
    const std::string message = context + problem.what();
    if (dynamic_cast<const InputError*>(&problem) != nullptr)
        Logger::instance().write_line(LogLevel::error, message);
    else
        Logger::instance().write(LogLevel::error, message);
}

template <typename... Args>
[[noreturn]] void fatal(Args&&... args)
{
    Logger::instance().write(LogLevel::fatal, std::forward<Args>(args)...);
    std::abort();
}

} // namespace niceops
