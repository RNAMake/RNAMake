#ifndef PYBIND11_BASE_HPP
#define PYBIND11_BASE_HPP
// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/cast.h>
#include <pybind11/operators.h>
#include <memory>

// base includes
#include <base/application.hpp>
#include <base/backtrace.h>
#include <base/cl_option.h>
#include <base/command_line_parser.hpp>
#include <base/env_manager.h>
#include <base/exception.h>
#include <base/file_io.h>
#include <base/log.h>
#include <base/option.h>
#include <base/settings.h>
#include <base/string.h>
#include <base/sys_interface.h>
#include <base/vector_container.h>


namespace base {
    // trampoline class

    class PyApplication : public Application {
    public:
        PyApplication() : Application() {
        }
    public:
        virtual
        void
        setup_options() override {
            throw ApplicationException("RNAMake.base.setup_options() must be overridden");
        }

    public:
        virtual
        void
        run() override {
            throw ApplicationException("RNAMake.base.run() must be overridden");
        }

    public:
        CommandLineOptions const&
        cl_options() const {
            return cl_options_;
        }

    public:
        void
        cl_options(CommandLineOptions & cl_opts) {
            cl_options_ = cl_opts;
        }

    public:
        CommandLineParser const&
        cl_parser() const {
            return cl_parser_;
        }

    public:
        void
        cl_parser(CommandLineParser & cl_parser ) {
            cl_parser_ = cl_parser;
        }
    };



    namespace py = pybind11;

    void
    add_bindings(pybind11::module_ &m) {
        // Exceptions
        py::register_exception<RNAMakeException>(m, "RNAMakeException");
        py::register_exception<RNAMakeIOException>(m, "RNAMakeIOException");
        py::register_exception<RNAMakeImplementationExcepetion>(m, "RNAMakeImplementationExcepetion");
        py::register_exception<CommandLineOptionException>(m, "CommandLineOptionException");
        py::register_exception<ApplicationException>(m, "ApplicationException");
        py::register_exception<OptionException>(m, "OptionException");

        // free functions
        m.def("base_dir", [](String const &path) -> String {
                  return base_dir(path);
              },
              py::arg("path")
        );

        m.def("demangle", [](std::string &str) -> std::string {
                  return demangle(str);
              },
              py::arg("str")
        );

        m.def("determine_string_data_type", [](String const &str) -> DataType {
                  return determine_string_data_type(str);
              },
              py::arg("str")
        );

        m.def("execute_command", [](const char *cmd) -> String {
                  return execute_command(cmd);
              },
              py::arg("cmd")
        );

        m.def("execute_command_json", [](const char *cmd) -> nlohmann::json {
                  return execute_command_json(cmd);
              },
              py::arg("cmd")
        );

        m.def("file_exists", [](String const &name) -> bool {
                  return file_exists(name);
              },
              py::arg("name")
        );

        m.def("filename", [](String const &path) -> String {
                  return filename(path);
              },
              py::arg("path")
        );

        m.def("get_lines_from_file", [](String const &name) -> Strings {
                  return get_lines_from_file(name);
              },
              py::arg("name")
        );

        m.def("get_os_name", []() -> String {
                  return get_os_name();
              }
        );

        py::enum_<LogLevel>(m, "LogLevel")
                .value("DEBUG", LogLevel::DEBUG)
                .value("ERROR", LogLevel::ERROR)
                .value("FATAL", LogLevel::FATAL)
                .value("INFO", LogLevel::INFO)
                .value("VERBOSE", LogLevel::VERBOSE)
                .value("WARN", LogLevel::WARN);

        m.def("init_logging", [](LogLevel log_level) { // TODO add the log level stuff
                  init_logging(log_level);
              },
              py::arg("log_level") = LogLevel::INFO
        );

        m.def("is_dir", [](String const &path) -> int {
                  return is_dir(path);
              },
              py::arg("path")
        );

        m.def("is_number", [](String const &s) -> bool {
                  return is_number(s);
              },
              py::arg("s")
        );

        m.def("join_by_delimiter", [](Strings const &strs, String const &delimiter) -> String {
                  return join_by_delimiter(strs, delimiter);
              },
              py::arg("strs"),
              py::arg("delimiter")
        );

        m.def("lib_path", []() -> String {
                  return lib_path();
              }
        );

        m.def("log_level_from_str", [](String const &s) -> LogLevel {
                  return log_level_from_str(s);
              },
              py::arg("s")
        );

        m.def("ltrim", [](String &s) -> String & {
                  return ltrim(s);
              },
              py::arg("s")
        );

        m.def("motif_dirs", []() -> String {
                  return motif_dirs();
              }
        );

        m.def("print_backtrace", []() {
                  print_backtrace();
              }
        );

        m.def("replace_all", [](String &context, String const &from, String const &to) -> String & {
                  return replace_all(context, from, to);
              },
              py::arg("context"),
              py::arg("from"),
              py::arg("to")
        );

        m.def("resources_path", []() -> String {
                  return resources_path();
              }
        );

        m.def("rtrim", [](String &s) -> String & {
                  return rtrim(s);
              },
              py::arg("s")
        );

        m.def("save_backtrace", []() {
                  save_backtrace();
              }
        );

        m.def("split_str_by_delimiter", [](String& s, String& delimiter) -> Strings {
                  return split_str_by_delimiter(s, delimiter);
              },
              py::arg("s"),
              py::arg("delimiter")
        );

        m.def("tokenize_line", [](String const &raw_line) -> Strings {
                  return tokenize_line(raw_line);
              },
              py::arg("raw_line")
        );

        m.def("trim", [](String &s) -> String & {
                  return trim(s);
              },
              py::arg("s")
        );

        m.def("unittest_resource_dir", []() -> String {
                  return unittest_resource_dir();
              }
        );

        m.def("x3dna_path", []() -> String {
                  return x3dna_path();
              }
        );
        // classes
        py::class_<PyApplication, std::shared_ptr<PyApplication>>(m, "Application")
		        // ctors
		        .def(py::init<>())
		        // methods
		        .def("setup_options",[] (PyApplication  & ptr) { ptr.setup_options(); } )
		        .def("parse_command_line",[] (PyApplication  & ptr, int argc, Strings const& vec) {
                    std::vector<char *> ptrs; ptrs.reserve(vec.size()); for (auto &v : vec) { ptrs.push_back(const_cast<char *>(v.c_str())); };
                    ptr.parse_command_line(argc, (const char **) (&ptrs[0])); },
                    py::arg("argc"),
                    py::arg("argv")
                )
		        .def("run",[] (PyApplication  & ptr) { ptr.run(); } )
		        .def("add_option",[] (PyApplication  & ptr, String const & name, int const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); },
                    py::arg("name"),
                    py::arg("val"),
                    py::arg("type"),
                    py::arg("required") = false
                )
                .def("add_option",[] (PyApplication  & ptr, String const & name, float const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); },
                    py::arg("name"),
                    py::arg("val"),
                    py::arg("type"),
                    py::arg("required") = false
                )
                .def("add_option",[] (PyApplication  & ptr, String const & name, bool const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); },
                    py::arg("name"),
                    py::arg("val"),
                    py::arg("type"),
                    py::arg("required") = false
                )
                .def("add_option",[] (PyApplication  & ptr, String const & name, String const & val, OptionType const & type, bool required) { ptr.add_option(name, val, type, required); },
                    py::arg("name"),
                    py::arg("val"),
                    py::arg("type"),
                    py::arg("required") = false)
		        .def("add_cl_options",[] (PyApplication  & ptr, Options const & opts, String prefix) { ptr.add_cl_options(opts, prefix); }, py::arg("opts"), py::arg("prefix") )
		        .def("get_int_option",[] (PyApplication  & ptr, String const & name) ->  int { return ptr.get_int_option(name); }, py::arg("name") )
		        .def("get_float_option",[] (PyApplication  & ptr, String const & name) ->  float { return ptr.get_float_option(name); } , py::arg("name"))
		        .def("get_string_option",[] (PyApplication  & ptr, String const & name) ->  String { return ptr.get_string_option(name); } , py::arg("name"))
		        .def("get_bool_option",[] (PyApplication  & ptr, String const & name) ->  bool { return ptr.get_bool_option(name); } , py::arg("name"))
	            .def("cl_options", [] (PyApplication & ptr) { return ptr.cl_options(); })
                .def("cl_options", [] (PyApplication & ptr, CommandLineOptions & clopts) { ptr.cl_options(clopts); }, py::arg("clopts"))
                .def("cl_parser",  [] (PyApplication & ptr) { return ptr.cl_parser(); })
                .def("cl_parser",  [] (PyApplication & ptr, CommandLineParser & clparser) { ptr.cl_parser(clparser); }, py::arg("clparser"))
		        ;

        py::class_<CommandLineOption, std::shared_ptr<CommandLineOption>>(m, "CommandLineOption")
                // ctors
                .def(py::init<String const &, int const &, OptionType const &, bool>(), py::arg("name"), py::arg("value"), py::arg("type"), py::arg("required"))
                .def(py::init<String const &, float const &, OptionType const &, bool>(), py::arg("name"), py::arg("value"), py::arg("type"), py::arg("required"))
                .def(py::init<String const &, bool const &, OptionType const &, bool>(), py::arg("name"), py::arg("value"), py::arg("type"), py::arg("required"))
                .def(py::init<String const &, String const &, OptionType const &, bool>(), py::arg("name"), py::arg("value"), py::arg("type"), py::arg("required"))
                .def(py::init<Option const &>(), py::arg("option"))
                // methods
                .def("filled", [](CommandLineOption &ptr) -> bool { return ptr.filled(); })
                .def("required", [](CommandLineOption &ptr) -> bool { return ptr.required(); })
                .def("filled", [](CommandLineOption &ptr, bool const &filled) { ptr.filled(filled); })
                // inherited methods
                .def("get_float", [](Option &ptr) -> float { return ptr.get_float(); })
                .def("get_int", [](Option &ptr) -> int { return ptr.get_int(); })
                .def("get_string", [](Option &ptr) -> String const & { return ptr.get_string(); })
                .def("get_bool", [](Option &ptr) -> bool const & { return ptr.get_bool(); })
                .def("value", [](Option &ptr, int const &i) { ptr.value(i); }, py::arg("i"))
                .def("value", [](Option &ptr, char const *c_str) { ptr.value(c_str); }, py::arg("c_str"))
                .def("value", [](Option &ptr, bool const &boolean) { ptr.value(boolean); }, py::arg("boolean"))
                .def("value", [](Option &ptr, float const &real) { ptr.value(real); }, py::arg("real"))
                .def("value", [](Option &ptr, String const &str) { ptr.value(str); }, py::arg("str"))
                .def("name", [](Option const &ptr) -> String const & { return ptr.name(); })
                .def("type", [](Option const &ptr) -> OptionType const & { return ptr.type(); })
                .def("type_name", [](Option &ptr) -> String { return ptr.type_name(); });

        py::class_<CommandLineOptions, std::shared_ptr<CommandLineOptions>>(m, "CommandLineOptions")
                // ctors
                .def(py::init<>())
                        // methods
                .def("begin", [](CommandLineOptions &ptr) -> std::vector<CommandLineOptionOP>::iterator {
                    return ptr.begin();
                })
                .def("end", [](CommandLineOptions &ptr) -> std::vector<CommandLineOptionOP>::iterator {
                    return ptr.end();
                })
                .def("begin", [](CommandLineOptions const &ptr) -> std::vector<CommandLineOptionOP>::const_iterator {
                    return ptr.begin();
                })
                .def("end", [](CommandLineOptions const &ptr) -> std::vector<CommandLineOptionOP>::const_iterator {
                    return ptr.end();
                })
                .def("add_option",
                     [](CommandLineOptions &ptr, String const &name, String const &value, OptionType const &type,
                        bool required) { ptr.add_option(name, value, type, required); },
                        py::arg("name"),
                        py::arg("value"),
                        py::arg("type"),
                        py::arg("required") = false)
                .def("add_option",
                     [](CommandLineOptions &ptr, String const &name, int const &value, OptionType const &type,
                        bool required) { ptr.add_option(name, value, type, required); },
                        py::arg("name"),
                        py::arg("value"),
                        py::arg("type"),
                        py::arg("required") = false)
                .def("add_option",
                     [](CommandLineOptions &ptr, String const &name, bool const &value, OptionType const &type,
                        bool required) { ptr.add_option(name, value, type, required); },
                        py::arg("name"),
                        py::arg("value"),
                        py::arg("type"),
                        py::arg("required") = false)
                .def("add_option",
                     [](CommandLineOptions &ptr, String const &name, float const &value, OptionType const &type,
                        bool required) { ptr.add_option(name, value, type, required); },
                        py::arg("name"),
                        py::arg("value"),
                        py::arg("type"),
                        py::arg("required") = false)
                .def("add_options", [](CommandLineOptions &ptr, Options const &opts) {
                    ptr.add_options(opts);
                }, py::arg("opts"))
                .def("parse_command_line", [](CommandLineOptions &ptr, const int argc, Strings &vec) {
                    std::vector<char *> ptrs;
                    ptrs.reserve(vec.size());
                    for (auto &v : vec) { ptrs.push_back(const_cast<char *>(v.c_str())); };
                    ptr.parse_command_line(argc, (const char **) (&ptrs[0]));
                },py::arg("argc"), py::arg("argv"))
                .def("get_int", [](CommandLineOptions const &ptr, String const &name) -> float { return ptr.get_int(name); })
                .def("get_float", [](CommandLineOptions const &ptr, String const &name) -> float { return ptr.get_float(name); })
                .def("get_string", [](CommandLineOptions const &ptr, String const &name) -> String { return ptr.get_string(name); })
                .def("get_bool", [](CommandLineOptions const &ptr, String const &name) -> bool {return ptr.get_bool(name); })
                .def("has_option", [](CommandLineOptions const &ptr, String const &name) -> bool { return ptr.has_option(name); })
                .def("set_value",[](CommandLineOptions &ptr, String const &name, String const &val) { ptr.set_value(name, val); }, py::arg("name"), py::arg("value"))
                .def("set_value", [](CommandLineOptions &ptr, String const &name, int const &val) { ptr.set_value(name, val); }, py::arg("name"), py::arg("value"))
                .def("set_value", [](CommandLineOptions &ptr, String const &name, bool const &val) { ptr.set_value(name, val); }, py::arg("name"), py::arg("value"))
                .def("set_value", [](CommandLineOptions &ptr, String const &name, float const &val) { ptr.set_value(name, val); }, py::arg("name"), py::arg("value"))
                .def("is_filled", [](CommandLineOptions const &ptr, String const &name) -> bool {return ptr.is_filled(name);}, py::arg("name"));

        py::class_<CommandLineParser, std::shared_ptr<CommandLineParser>>(m, "CommandLineParser")
                // ctors
                .def(py::init<>())
                        // methods
                .def("assign_options",
                     [](CommandLineParser &ptr, CommandLineOptions const &cl_options, Options &options, String& prefix) {
                         ptr.assign_options(cl_options, options, prefix);
                     },py::arg("cl_options"), py::arg("options"), py::arg("prefix") = "");

        py::class_<EnvManager, std::shared_ptr<EnvManager>>(m, "EnvManager")
                // ctors
                .def(py::init<Strings const &>(), py::arg("env_vars"))
                        // methods
                .def("add_env", [](EnvManager &ptr, String const &env) {
                    ptr.add_env(env);
                },py::arg("env"))
                .def("set_envs", [](EnvManager &ptr) {
                    ptr.set_envs();
                });


        py::class_<Option, std::shared_ptr<Option>>(m, "Option")
                // ctors
                .def(py::init<>())
                .def(py::init<Option const &>(), py::arg("opt"))
                .def(py::init<String const &, float const &, OptionType const &>(),
                        py::arg("name"), py::arg("value"), py::arg("type"))
                .def(py::init<String const &, String const &, OptionType const &>(),
                     py::arg("name"), py::arg("value"), py::arg("type"))
                .def(py::init<String const &, char const *, OptionType const &>(),
                     py::arg("name"), py::arg("value"), py::arg("type"))
                .def(py::init<char const *, char const *, OptionType const &>(),
                     py::arg("name"), py::arg("value"), py::arg("type"))
                .def(py::init<String const &, bool const &, OptionType const &>(),
                py::arg("name"), py::arg("value"), py::arg("type"))
                .def(py::init<String const &, int const &, OptionType const &>(),
                py::arg("name"), py::arg("value"), py::arg("type"))
                        // methods
                .def("get_float", [](Option &ptr) -> float {
                    return ptr.get_float();
                })
                .def("get_int", [](Option &ptr) -> int {
                    return ptr.get_int();
                })
                .def("get_string", [](Option &ptr) -> String const & {
                    return ptr.get_string();
                })
                .def("get_bool", [](Option &ptr) -> bool const & {
                    return ptr.get_bool();
                })
                .def("value", [](Option &ptr, int const &i) {
                    ptr.value(i);
                }, py::arg("i"))
                .def("value", [](Option &ptr, char const *c_str) {
                    ptr.value(c_str);
                }, py::arg("c_str"))
                .def("value", [](Option &ptr, int const &i) {
                    ptr.value(i);
                }, py::arg("i"))
                .def("value", [](Option &ptr, int const &i) {
                    ptr.value(i);
                }, py::arg("i"))
                .def("value", [](Option &ptr, int const &i) {
                    ptr.value(i);
                }, py::arg("i"))
                .def("name", [](Option const &ptr) -> String const & {
                    return ptr.name();
                })
                .def("type", [](Option const &ptr) -> OptionType const & {
                    return ptr.type();
                })
                .def("type_name", [](Option &ptr) -> String {
                    return ptr.type_name();
                });

        py::class_<OptionClass, std::shared_ptr<OptionClass>>(m, "OptionClass");

        py::enum_<OptionType>(m, "OptionType")
            .value("BOOL", OptionType::BOOL)
            .value("INT", OptionType::INT)
            .value("STRING", OptionType::STRING)
            .value("FLOAT",OptionType::FLOAT)
            ;

        py::class_<Options, std::shared_ptr<Options>>(m, "Options")
                // ctors
                .def(py::init<>())
                .def(py::init<String const &>(), py::arg("name"))
                .def(py::init<Options const &>(), py::arg("opts"))
                        // methods
                .def("begin", [](Options &ptr) -> std::vector<OptionOP>::iterator {
                    return ptr.begin();
                })
                .def("end", [](Options &ptr) -> std::vector<OptionOP>::iterator {
                    return ptr.end();
                })
                .def("begin", [](Options const &ptr) -> std::vector<OptionOP>::const_iterator {
                    return ptr.begin();
                })
                .def("end", [](Options const &ptr) -> std::vector<OptionOP>::const_iterator {
                    return ptr.end();
                })
                .def("size", [](Options &ptr) -> size_t {
                    return ptr.size();
                })
                .def("lock_option_adding", [](Options &ptr) {
                    ptr.lock_option_adding();
                })
                .def("add_option", [](Options &ptr, String const &name, int const &val, OptionType const &type) {
                    ptr.add_option(name, val, type);
                })
                .def("add_option", [](Options &ptr, String const &name, bool const &val, OptionType const &type) {
                    ptr.add_option(name, val, type);
                })
                .def("add_option", [](Options &ptr, String const &name, String const &val, OptionType const &type) {
                    ptr.add_option(name, val, type);
                })
                .def("add_option", [](Options &ptr, String const &name, float const &val, OptionType const &type) {
                    ptr.add_option(name, val, type);
                })
                .def("get_int", [](Options &ptr, String const &name) -> int {
                    return ptr.get_int(name);
                })
                .def("get_float", [](Options &ptr, String const &name) -> float {
                    return ptr.get_float(name);
                })
                .def("get_string", [](Options &ptr, String const &name) -> String const & {
                    return ptr.get_string(name);
                })
                .def("get_bool", [](Options &ptr, String const &name) -> bool {
                    return ptr.get_bool(name);
                })
                .def("has_option", [](Options &ptr, String const &name) -> bool {
                    return ptr.has_option(name);
                })
                .def("set_value", [](Options &ptr, String const &name, bool const &val) { ptr.set_value(name, val); })
                .def("set_value", [](Options &ptr, String const &name, float const &val) { ptr.set_value(name, val); })
                .def("set_value", [](Options &ptr, String const &name, int const &val) { ptr.set_value(name, val); })
                .def("set_value",
                     [](Options &ptr, String const &name, String const &val) { ptr.set_value(name, val); });


/*
        py::class_<VectorContainer, std::shared_ptr<VectorContainer>>(m, "VectorContainer")
		// ctors
		.def(py::init<std::vector<T> const &>())
		// methods
		.def("size",[] (VectorContainer  & ptr) ->  size_t {
		 return ptr.size(); } )
		.def("begin",[] (VectorContainer const & ptr) -> const_iterator {
		 return ptr.begin(); } )
		.def("end",[] (VectorContainer const & ptr) -> const_iterator {
		 return ptr.end(); } )
		.def("at",[] (VectorContainer const & ptr, Index i) ->  T const & {
		 return ptr.at(i); } )
		.def("get_data",[] (VectorContainer  & ptr) ->  std::vector<T> const & {
		 return ptr.get_data(); } )
		// operators
		.def(py::self [] Index)
		;
*/
    }
}
#endif // PYBIND11_BASE_HPP