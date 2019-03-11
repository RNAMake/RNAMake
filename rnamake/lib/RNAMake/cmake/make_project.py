import os
import fnmatch
import argparse

depends = {
    'base' : [],
    'math' : ['base'],
    'data_structure' : ['base'],
    'util' : ['math'],
    'vienna' : ['base'],
    'secondary_structure' : ['util'],
    'eternabot' : ['vienna', 'secondary_structure'],
    'structure' : ['util'],
    'motif' : ['structure', 'secondary_structure'],
    'motif_tools' : ['motif'],
    'resources' : ['motif'],
    'motif_data_structures' : ['resources', 'data_structure'],
    'thermo_fluctuation' : ['motif_data_structures'],
    'motif_state_search' : ['motif_data_structures'],
    'sequence_optimizer' : ['motif_data_structures', 'eternabot'],
    'all' : ['motif_tools', 'thermo_fluctuation', 'motif_state_search',
                  'sequence_optimizer']
}

libs = "base math data_structure util vienna secondary_structure eternabot structure motif motif_tools resources motif_data_structures thermo_fluctuation motif_state_search sequence_optimizer"
#libs = "base math data_structure util vienna secondary_structure eternabot structure motif motif_tools "
all_lib_paths = libs.split()


file_path = os.path.realpath(__file__)
file_path = os.path.realpath(__file__)
spl = file_path.split("/")
base_dir = "/".join(spl[:-2])

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-no_unittests', action="store_true")
    parser.add_argument('-no_apps', default="store_false", action="store_true")
    parser.add_argument('-last_ta', action="store_true")

    args = parser.parse_args()

    return args

def get_lib_paths(args):
    return all_lib_paths

def get_cmake_lists_header():
    s =  "cmake_minimum_required(VERSION 2.8.12)\n"
    s += "project(rnamake_new)\n\n"
    s += "set(CMAKE_BUILD_TYPE Release)\n"
    s += "include(%s)\n\n" % (base_dir + "/cmake/build/compiler.cmake")
    s += "# Include path for Python header files\n"
    s += "# Include paths for RNAMake src\n"
    s += "include_directories(%s)\n" % (base_dir + "/src/")
    s += "include_directories(%s)\n" % (base_dir + "/src/plog/")
    s += "include_directories(%s)\n\n" % (base_dir + "/unittests/")
    s += "include_directories(%s)\n\n" % (base_dir + "/apps/")
    s += "# sqlite libraries\n"
    s += "find_library(SQLITE3_LIBRARY sqlite3)\n\n"
    return s

def get_pretty_lib_name(name):
    s = ""
    for i in range(75):
        s += "#"
    s += "\n"
    s += "# " + name + "\n"
    for i in range(75):
        s += "#"
    s += "\n\n"
    return s

def get_lib_file_declaration(lib):
    abs_path = base_dir + "/src/" + lib + "/"
    cc_files = get_cc_files_in_dir(abs_path)
    s  = "set(%s_files\n" % (lib)
    for cc_file in cc_files:
        s += "\t" + cc_file + "\n"
    s += ")\n\n"
    return s

def get_unittests_apps_for_library(lib):
    unittest_apps = []
    abs_path = base_dir + "/unittests/" + lib + "/"

    if not os.path.isdir(abs_path):
        return ""

    unittest_files = get_cc_files_in_dir(abs_path)

    s = ""
    for unit in unittest_files:
        spl = unit.split("/")[-1].split(".")
        prog_name = spl[0]

        s += "add_executable(" + prog_name + " " +  unit + ")\n"
        s += "target_link_libraries(" + prog_name + " %s_lib)\n\n" % lib

    return s

def get_build_library_declaration(lib):
    s  = "add_library(%s ${%s})\n" % (lib + "_lib", lib + "_files")
    s += get_linking_declaration(lib)
    return s

def get_linking_declaration(lib):
    s = "target_link_libraries(%s" % (lib + "_lib")
    for depend in depends[lib]:
        s += " " + depend + "_lib "
    if lib == "util":
        s += " ${SQLITE3_LIBRARY})\n\n"
    else:
        s += ")\n\n"
    return s

def get_cc_files_in_dir(path):
    cpp_files = []
    for root, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename[-2:] == 'cc' or filename[-3:] == 'cpp':
                    final_root = root
                    if final_root[-1] != "/":
                        final_root += "/"
                    cpp_files.append(final_root + filename)
    return cpp_files

def get_applications():
    path = base_dir + "/cmake/build/apps.txt"
    f = open(path)
    lines = f.readlines()
    f.close()

    s = ""
    for l in lines:
        if l[0] == "#":
         continue

        spl = l.split()
        app_dir = base_dir + "/apps/" + spl[0]
        cpp_files = get_cc_files_in_dir(app_dir)

        symlink = spl[0]+"_symlink"
        s += "add_executable("+spl[0] + " "
        for f in cpp_files:
            s += f + " "
        s += ")\n"
        s += "target_link_libraries("+spl[0]+" all_lib)\n"
        s += "\n\n"
    return s

def get_unittests_for_dir(d_name):
    unittest_apps = []

    if not os.path.isdir('../../unittests/'+d_name):
        return []

    for root, dirnames, filenames in os.walk('../../unittests/'+d_name):
        for filename in fnmatch.filter(filenames, '*.c*'):
            path = os.path.join(root, filename)

            f = open(path)
            fail = 0
            for l in f.readlines():
                if len(l) < 8:
                    continue
                if l[:9] == 'TEST_CASE':
                    unittest_apps.append((os.path.join(root, filename)))
                    fail = 1
                    break

    return unittest_apps

def write_cmake_lists(path, args):
    f = open(path+"/CMakeLists.txt", "w")
    f.write(get_cmake_lists_header())

    lib_paths = get_lib_paths(args)
    for lib in lib_paths:
        f.write(get_pretty_lib_name(lib))
        f.write(get_lib_file_declaration(lib))
        f.write(get_build_library_declaration(lib))
        if not args.no_unittests:
            f.write(get_unittests_apps_for_library(lib))

    #f.write(get_pretty_lib_name("all"))
    #f.write("add_library(all_lib " + base_dir + "/src/main.cpp)\n")
    #f.write(get_linking_declaration("all"))

    if not args.no_apps:
        f.write(get_applications())
    f.close()

if __name__ == '__main__':
    args = parse_args()
    write_cmake_lists("./", args)














