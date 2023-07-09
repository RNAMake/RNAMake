def is_header_file(file_name):
    """Method that checks if the file is a C++ header file"""
    file_ending = file_name.split('.')[-1].lower()
    for ending in "h|hh|hpp|hxx|h++".split('|'):
        if ending == file_ending:
            return True
    return False

def is_source_file(file_name):
    """Method that checks if the file is a C++ source file"""
    file_ending = file_name.split('.')[-1].lower()
    for ending in "c|cc|cpp|cxx|c++".split('|'):
        if ending == file_ending:
            return True
    return False

def is_cpp_file(file_name):
    """Method that checks if the file is a C++ file"""
    return is_header_file(file_name) or is_source_file(file_name)

def get_file_name(raw_name):
    """Method that gets the file_name from te input raw_name"""
    return raw_name.split('.')[-2]

def make_file_list(raw_file_names):
    """Method that takes a list of files and returns the singular list of files for the CMakeLists.txt file"""
    source_files = dict()
    header_files = dict()
    # looping through the file names and adding them to the set
    for file_name in raw_file_names:
        # need to check if the file is a C++ file
        if is_cpp_file(file_name):
            # sort into the appropriate dict
            name = get_file_name(file_name)
            if is_source_file(file_name):
                source_files[name] = file_name
            else:
                header_files[name] = file_name
    # all source files go into the output file set
    output_file_set = set(source_files.values())
    # and the header files that are not already included as source
    for name, file_path in header_files.items():
        if name not in source_files:
            output_file_set.add(file_path)
    return list(output_file_set)



if __name__ == "__main__":
    pass
   # file_names = [
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/file_io.cc",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/file_io.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/string.cc",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/command_line_parser.hpp",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/settings.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/backtrace.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/types.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/application.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/settings.cc",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/cl_option.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/option.h",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/application.cpp",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/util.hpp",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/cl_option.cc",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/command_line_parser.cpp",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/backtrace.cpp",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/option.cc",
   #     "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/string.h"
   # ]
   # output = make_file_list(file_names)

   # targets = {
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/file_io.cc",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/string.cc",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/types.h",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/settings.cc",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/application.cpp",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/util.hpp",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/cl_option.cc",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/command_line_parser.cpp",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/backtrace.cpp",
   #         "/Users/cjurich/projects/RNAMake/rnamake/lib/RNAMake/src/base/option.cc"
   # }

   # for f in output:
   #     print(f in targets)

   # print(len(set(output)),len(targets))h
