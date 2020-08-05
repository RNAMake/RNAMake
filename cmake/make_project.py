import os
import Utils
import glob
import re
from pathlib import Path
import argparse

def build_libraries(lib_list, dependencies, base_dir,static):
    """Method that builds out the libraries, taking a list of libraries and a dependency dictionary"""
    library_declarations = "#"*100 + "\n# Library declarations \n" + "#"*100 + '\n'
    for library in lib_list:
        library_declarations += "#"*75 + "\n# {NAME}\n".format(NAME=library) + "#"*75 + '\n'
        library_declarations += "set({LIB}_files\n\t{FILES} \n)\n".format(LIB=library, FILES="\n\t".join(Utils.make_file_list(
            [str(fp) for fp in Path(base_dir + "/" + library ).rglob("*")])
        ))
        #library_declarations += "add_library({LIB}_lib STATIC ${{{LIB}_files}})\n".format(LIB=library)
        library_declarations += "add_library({LIB}_lib {BUILD} ${{{LIB}_files}})\n".format(
                LIB=library,
                BUILD="STATIC" if static else ""
                )
        #library_declarations += "target_link_libraries({LIB}_lib {DEPENDS} -static)\n".format(
        library_declarations += "target_link_libraries({LIB}_lib {DEPENDS} {BUILD})\n".format(
            LIB=library,
            DEPENDS=' '.join([lib + "_lib" for lib in dependencies[library]] + ["${SQLITE3_LIBRARY} \"-ldl\"" if library == "util" else ""]),
            BUILD="-static" if static else ""
        )
    #library_declarations += "add_library(all_lib STATIC {DIR}/main.cpp)\ntarget_link_libraries(all_lib {LIBS} -static)\n".format(
    library_declarations += "add_library(all_lib {BUILD} {DIR}/main.cpp)\ntarget_link_libraries(all_lib {LIBS} {LINK})\n".format(
        LIBS=' '.join([lib+"_lib" for lib in dependencies["all"]]),DIR=base_dir,
        BUILD="STATIC" if static else "",
        LINK="-static" if static else ""
    )
    library_declarations += "#"*100 + '\n'
    return library_declarations

def build_unittests(lib_list, base_dir,static):
    """Method that builds out the declarations for the unittests"""
    unittest_declarations = "#"*100 + "\n# Unittest Declarations \n" + "#"*100 + '\n'
    for library in lib_list:
        unittest_declarations += "#"*75 + "\n# {MODULE} Tests \n".format(MODULE=library) + "#"*75 + '\n'
        for source_file in Utils.make_file_list(glob.glob(base_dir + "/" + library + "/*")):
            unittest_name = source_file.split('/')[-1].split('.')[0]
            unittest_declarations += "\tadd_executable({TEST} {SRC})\n".format(
                TEST=unittest_name, SRC=source_file
            )
            #unittest_declarations += "\ttarget_link_libraries({TEST} {LIB}_lib -static)\n".format(
            unittest_declarations += "\ttarget_link_libraries({TEST} {LIB}_lib {BUILD} )\n".format(
                TEST=unittest_name,
                LIB=library,
                BUILD="-static" if static else ""

            )
    unittest_declarations += "#"*100 + '\n'
    return unittest_declarations

def write_CML_file(content_list):
    """Method that exports the contents to the CMakeLists.txt file"""
    with open("CMakeLists.txt", "w")  as outfile:
        for content in content_list:
            outfile.write(content)

def build_apps(base_dir,static):
    """Method that creates declarations for the applications"""
    application_text = "#"*100 + "\n" + "# App declarations\n" + "#"*100
    with open("apps.txt", "r") as infile:
        for app_declaration in infile.readlines():
            if app_declaration.find("#") == -1:
                app_tokens = app_declaration.split()
                source_file = re.sub("\.\./","",app_tokens[-1])
                application_text += "#"*75 + "\n# {NAME}".format(NAME=app_tokens[0]) + "\n" + "#"*75 + "\n"
                application_text += "\tadd_executable({NAME} {SRC})\n".format(NAME=app_tokens[0],SRC="/".join([base_dir,source_file]))
                application_text += "\ttarget_link_libraries({NAME} all_lib {LINK} )\n".format(
                        NAME=app_tokens[0],
                        LINK="-static" if static else ""
                        )
                #application_text += "\ttarget_link_libraries({NAME} all_lib -static)\n".format(NAME=app_tokens[0])
    application_text += "#"*100 + '\n'
    return application_text

def build_header(base_dir,static):
    """Method that creates the header for the CMakeLists.txt file"""
    header_contents = "#"*100 + "\n# Project Level Info\n" + "#"*100 + '\n'
    header_contents +=  "cmake_minimum_required(VERSION 3.0)\n"
    header_contents+= "set(CMAKE_BUILD_TYPE Release)\n"
    header_contents+= "project(rnamake_new)\n\n"
    header_contents += "set(CMAKE_CXX_STANDARD 14)\n"
    # header_contents += "set(CMAKE_C_COMPILER clang)\n"
    # header_contents += "set(CMAKE_CXX_COMPILER clang++)\n"
    header_contents+= "set( CMAKE_CXX_FLAGS \" -pthread -L/opt/local/lib \" )\n"
    if static == True:
        header_contents+= "set(CMAKE_SHARED_LINKER_FLAGS \"-Wl,--no-as-needed -ldl\")\n"
    if static == True:
        header_contents+= "set(CMAKE_EXE_LINKER_FLAGS \" -lstdc++ -Wl,--no-as-needed -ldl \")\n"
    # header_contents+= "set(CMAKE_HOST_SYSTEM_VERSION 2.5)\n"
    header_contents+= "include({CMAKE})\n\n".format(CMAKE=(base_dir + "/cmake/build/compiler.cmake"))
    header_contents+= "# Include path for Python header files\n"
    header_contents+= "# Include paths for RNAMake src\n"
    header_contents+= "include_directories({DIR})\n".format(DIR=(base_dir + "/src/"))
    header_contents+= "include_directories({DIR})\n".format(DIR=base_dir + "/src/plog/")
    header_contents+= "include_directories({DIR})\n\n".format(DIR=base_dir + "/unittests/")
    header_contents+= "include_directories({DIR})\n\n".format(DIR=base_dir + "/apps/")
    header_contents+= "# sqlite libraries\n"
    header_contents+= "find_library(SQLITE3_LIBRARY sqlite3 VARIANT {BUILD} )\n".format(
                BUILD="static" if static else ""
            )
    #header_contents+= "find_library(SQLITE3_LIBRARY sqlite3 VARIANT static)\n"
    header_contents += "#"*100 + '\n'
    return header_contents

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-target",type=str,required=True)

    args = parser.parse_args()

    if args.target == "mac":
        static=False
    elif args.target == "windows":
        pass
    elif args.target == "linux":
        static=True
    else:
        raise Exception("Invalid value for \"-target\" flag. Acceptable values are \"mac\",\"windows\" and \"linux\"")


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
    'motif_data_structure' : ['resources', 'data_structure'],
    'thermo_fluctuation' : ['motif_data_structure'],
    'motif_search' : ['motif_data_structure'],
    'sequence_optimization' : ['motif_data_structure', 'eternabot'],
    'all' : ['motif_tools', 'thermo_fluctuation', 'motif_search',
                  'sequence_optimization']
    }

    libs = "base math data_structure util vienna secondary_structure eternabot structure motif motif_tools resources motif_data_structure thermo_fluctuation motif_search sequence_optimization".split()

    write_CML_file(
        [
        build_header(
                os.getcwd().replace("/cmake/build",""),
                static
                ),
        build_libraries(
                libs,
                depends,os.getcwd().replace("/cmake/build","/src"),
                static
                ),
        build_unittests(
                libs,os.getcwd().replace("/cmake/build","/unittests"),
                static
                ),
        build_apps(
                os.getcwd().split("/cmake/build")[0],
                static
            )
        ]
        )