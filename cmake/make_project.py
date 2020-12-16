import os
import Utils
import glob
import re
from pathlib import Path
import argparse


def build_libraries(lib_list, dependencies, base_dir, static):
    """Method that builds out the libraries, taking a list of libraries and a dependency dictionary"""
    library_declarations = "#" * 100 + "\n# Library declarations \n" + "#" * 100 + "\n"
    for library in lib_list:
        library_declarations += (
            "#" * 75 + "\n# {NAME}\n".format(NAME=library) + "#" * 75 + "\n"
        )
        library_declarations += "set({LIB}_files\n\t{FILES} \n)\n".format(
            LIB=library,
            FILES="\n\t".join(
                Utils.make_file_list(
                    [str(fp) for fp in Path(base_dir + "/" + library).rglob("*")]
                )
            )
        )
        library_declarations += "add_library({LIB}_lib {BUILD} ${{{LIB}_files}})\n".format(
            LIB=library, BUILD="STATIC" if static else ""
        ) 
        library_declarations += "target_link_libraries({LIB}_lib {DEPENDS} {BUILD})\n".format(
            LIB=library,
            DEPENDS=" ".join(
                [lib + "_lib" for lib in dependencies[library]]
                + ["sqlite3" if library == "util" else ""]
            ),
            BUILD="-static" if static else "",
        )
    library_declarations += "add_library(all_lib {BUILD} {DIR}/main.cpp)\ntarget_link_libraries(all_lib {LIBS} {LINK})\n".format(
        LIBS=" ".join([lib + "_lib" for lib in dependencies["all"]]),
        DIR=base_dir,
        BUILD="STATIC" if static else "",
        LINK="-static" if static else "",
    )
    library_declarations += "#" * 100 + "\n"
    return library_declarations


def build_unittests(lib_list, base_dir, static):
    """Method that builds out the declarations for the unittests"""
    unittest_declarations = (
        "#" * 100 + "\n# Unittest Declarations \n" + "#" * 100 + "\n"
    )
    for library in lib_list:
        unittest_declarations += (
            "#" * 75 + "\n# {MODULE} Tests \n".format(MODULE=library) + "#" * 75 + "\n"
        )
        for source_file in Utils.make_file_list(
            glob.glob(base_dir + "/" + library + "/*")
        ):
            unittest_name = source_file.split("/")[-1].split(".")[0]
            unittest_declarations += "\tadd_executable({TEST} {SRC})\n".format(
                TEST=unittest_name, SRC=source_file
            )
            unittest_declarations += "\ttarget_link_libraries({TEST} {LIB}_lib {BUILD} )\n".format(
                TEST=unittest_name, LIB=library, BUILD="-static" if static else ""
            )
            unittest_declarations += "\tadd_test({TEST} {TEST})\n".format(
                TEST=unittest_name
            )
    unittest_declarations += "#" * 100 + "\n"
    return unittest_declarations


def write_CML_file(content_list):
    """Method that exports the contents to the CMakeLists.txt file"""
    with open("CMakeLists.txt", "w") as outfile:
        for content in content_list:
            outfile.write(content)


def build_apps(base_dir, static):
    """Method that creates declarations for the applications"""
    hundred_dashes = "#"*100
    seventyfive_dashes = "#"*75
    
    application_text = f"{hundred_dashes}\n# App declarations\n{hundred_dashes}\n"
    
    with open(base_dir + "/cmake/build/apps.txt", "r") as infile:
        
        for app_declaration in infile.readlines():
            
            if len(re.findall("^ #",app_declaration)) > 0:
                continue
            
            app_tokens = app_declaration.split()
            for full_path in  Utils.make_file_list(
                [str(fp) for fp in Path(f"{base_dir}/apps/{app_tokens[0]}/").rglob("*")]
                                                    ):
                source_file = re.sub("\.\./", "", full_path)
                app_name = source_file.split('/')[-1].split('.')[0]
                 
                application_text += f"{seventyfive_dashes}\n# {app_name}\n{seventyfive_dashes}\n"
                application_text += f"\tadd_executable( {app_name} {full_path})\n"
                application_text += "\ttarget_link_libraries({NAME} all_lib {LINK}  )\n".format(
                    NAME=app_name, LINK="-static" if static else ""
                )
    
    application_text += f"{hundred_dashes}\n"
    
    return application_text


def build_header(base_dir, static, target):
    """Method that creates the header for the CMakeLists.txt file"""
    header_contents = "#" * 100 + "\n# Project Level Info\n" + "#" * 100 + "\n"
    header_contents += "cmake_minimum_required(VERSION 3.0)\n"
    header_contents += "set(CMAKE_BUILD_TYPE Release)\n"
    header_contents += "project(RNAMake)\n"
    header_contents += "enable_testing()\n\n"

    if static == True:
        # header_contents+= "set(CMAKE_SHARED_LINKER_FLAGS \"-Wl,--no-as-needed -ldl\")\n"
        header_contents += 'set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed ")\n'
    if static == True:
        # header_contents+= "set(CMAKE_EXE_LINKER_FLAGS \" -lstdc++ -Wl,--no-as-needed,--no-export-dynamic -ldl \")\n"
        header_contents += 'set(CMAKE_EXE_LINKER_FLAGS " -lstdc++ -Wl,--no-as-needed,--no-export-dynamic ")\n'
    # header_contents+= "set(CMAKE_HOST_SYSTEM_VERSION 2.5)\n"
    header_contents += """include({CMAKE})
set(RNAMAKE {RNAMAKE})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories({BASE})
include_directories({EXTERN})
include_directories({UNITTESTS})
include_directories({APPS})
include({SQLITE})
 """.format(
        RNAMAKE=base_dir,
        SQLITE=base_dir + "/cmake/build/sqlite.cmake",
        CMAKE=base_dir + "/cmake/build/compiler.cmake",
        EXTERN=base_dir + "/src/external/",
        BASE=base_dir + "/src/",
        UNITTESTS=base_dir + "/unittests/",
        APPS=base_dir + "/apps/",
    )
    header_contents += "#" * 100 + "\n"
    return header_contents


def get_base_dir():
    """Method that returns the base dir of RNAMake"""
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    return "/".join(spl[:-2])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-target",
        help='enter the target system you are compiling for. Valid choices are: "mac", "windows", or "linux"',
        type=str,
        required=True,
    )

    args = parser.parse_args()

    if args.target == "mac":
        static = False
    elif args.target == "windows":
        static = True
    elif args.target == "linux":
        static = True
        # static=True #TODO this needs to change
    else:
        raise Exception(
            'Invalid value for "-target" flag. Acceptable values are "mac","windows" and "linux"'
        )

    depends = {
        "base": [],
        "math": ["base"],
        "data_structure": ["base"],
        "util": ["math"],
        "vienna": ["base"],
        "secondary_structure": ["util"],
        "eternabot": ["vienna", "secondary_structure"],
        "structure": ["util"],
        "motif": ["structure", "secondary_structure"],
        "motif_tools": ["motif"],
        "resources": ["motif"],
        "motif_data_structure": ["resources", "data_structure"],
        "thermo_fluctuation": ["motif_data_structure"],
        "motif_search": ["motif_data_structure"],
        "sequence_optimization": ["motif_data_structure", "eternabot"],
        "all": [
            "motif_tools",
            "thermo_fluctuation",
            "motif_search",
            "sequence_optimization",
        ],
    }

    libs = "base math data_structure util vienna secondary_structure eternabot structure motif motif_tools resources motif_data_structure thermo_fluctuation motif_search sequence_optimization".split()
    base_dir = get_base_dir()
    
    write_CML_file(
        [
            build_header(
                base_dir, static, args.target if args.target == "windows" else "NONE"
            ),
            build_libraries(libs, depends, base_dir + "/src", static),
            build_unittests(libs, base_dir + "/unittests", static),
            build_apps(base_dir, static),
        ]
    )
#
