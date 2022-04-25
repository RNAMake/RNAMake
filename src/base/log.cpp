//
// Created by Joseph Yesselman on 3/10/19.
//

#include <cctype>
#include <cwctype>

// RNAMake Headers
#include <base/log.hpp>

namespace base {

void init_logging(LogLevel log_level) {
  static plog::ColorConsoleAppender<plog::CustomFormatter> consoleAppender;
  plog::init((plog::Severity)log_level, &consoleAppender);
}

void init_logging_with_file(LogLevel log_level) {
  static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("logs.csv",
                                                                    0, 1);
  static plog::ColorConsoleAppender<plog::CustomFormatter> consoleAppender;
  plog::init((plog::Severity)log_level, &fileAppender)
      .addAppender(&consoleAppender);
}

LogLevel log_level_from_str(const String & s) {
  auto lower_str = s;
  for (char& p : lower_str) {
    p = std::towlower(p);
  }

  if (lower_str == "fatal") {
    return LogLevel::FATAL;
  }
#if defined(_WIN32) || defined(_WIN64)
  else if (lower_str == "error") {
    return LogLevel::WIN_ERROR;
  }
#else
  else if (lower_str == "error") {
    return LogLevel::ERROR;
  }
#endif
  else if (lower_str == "warn") {
    return LogLevel::WARN;
  } else if (lower_str == "info") {
    return LogLevel::INFO;
  } else if (lower_str == "debug") {
    return LogLevel::DEBUG;
  } else if (lower_str == "verbose") {
    return LogLevel::VERBOSE;
  } else {
    throw std::runtime_error(lower_str + " is not a recongized log level");
  }
}

}  // namespace base
