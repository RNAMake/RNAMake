//
// Created by Joseph Yesselman on 3/10/19.
//

#ifndef RNAMAKE_NEW_LOG_H
#define RNAMAKE_NEW_LOG_H

#include <plog/Log.h>
#include <plog/Appenders/ColorConsoleAppender.h>

#include <base/types.h>

namespace plog  {

class CustomFormatter {
public:
    static util::nstring header() { return util::nstring(); }

    static util::nstring format(const Record& record) {
        tm t;
        util::localtime_s(&t, &record.getTime().time);

        util::nostringstream ss;
        ss << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_hour << PLOG_NSTR(":") << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_min << PLOG_NSTR(":") << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_sec <<  PLOG_NSTR(" ");
        ss << std::setfill(PLOG_NSTR(' ')) << std::setw(5) << std::left << severityToString(record.getSeverity()) << PLOG_NSTR(" ");
        ss << PLOG_NSTR("[") << record.getFunc() << PLOG_NSTR("@") << record.getLine() << PLOG_NSTR("] ");
        ss << record.getMessage() << PLOG_NSTR("\n");

        return ss.str();
    }
};
}

namespace base {

// keep values in project
enum class LogLevel {
    FATAL   = plog::fatal,
    ERROR   = plog::error,
    WARN    = plog::warning,
    INFO    = plog::info,
    DEBUG   = plog::debug,
    VERBOSE = plog::verbose
};

LogLevel
log_level_from_str(
        String const &);

void
init_logging(
        LogLevel log_level = LogLevel::INFO);


}

#endif //RNAMAKE_NEW_LOG_H
