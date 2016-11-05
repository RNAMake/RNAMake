import os
import fnmatch

from rnamake import util, settings

BIN_PATH = settings.LIB_PATH + "/lib/RNAMake/bin/"
CMAKE_PATH = settings.LIB_PATH + "/lib/RNAMake/cmake/build/"

def get_cc_files_in_dir(path):
    cpp_files = []
    for root, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename[-2:] == 'cc' or filename[-3:] == 'cpp':
                    cpp_files.append(filename)
    return cpp_files

def write_application_cmake_file():

    fsum = open(CMAKE_PATH+"apps.cmake", "w")

    f = open(CMAKE_PATH+"apps.txt")
    lines = f.readlines()
    f.close()
    for l in lines:
        if l[0] == "#":
         continue

        spl = l.split()
        app_dir = '../../apps/'+spl[0]+"/"
        cpp_files = get_cc_files_in_dir(app_dir)

        symlink = spl[0]+"_symlink"
        fsum.write("add_executable("+spl[0] + " ")
        for f in cpp_files:
            fsum.write(app_dir + "/" + f + " ")
        fsum.write(")\n")
        fsum.write("target_link_libraries("+spl[0]+" all_libs)\n")
        fsum.write("add_custom_target("+symlink+" ALL)\n")
        fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND cmake -E create_symlink "+BIN_PATH+spl[0]+" "+ spl[0]+")\n")
        fsum.write("\n\n")

    fsum.close()


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


def get_integration_tests_for_dir(d_name):
    intergration_apps = []

    if not os.path.isdir('../../unittests/integration/'+d_name):
        return []

    for root, dirnames, filenames in os.walk('../../unittests/integration/'+d_name):
        for filename in fnmatch.filter(filenames, '*.c*'):
            path = os.path.join(root, filename)

            f = open(path)
            fail = 0
            for l in f.readlines():
                if len(l) < 8:
                    continue
                if l[:9] == 'TEST_CASE':
                    intergration_apps.append((os.path.join(root, filename)))
                    fail = 1
                    break

    return intergration_apps


c_dir = settings.LIB_PATH + "/lib/RNAMake/"

# make sure these directories exist
if not os.path.isdir(c_dir + "bin"):
    os.mkdir(c_dir + "bin")

if not os.path.isdir(c_dir + "bin/unittests"):
    os.mkdir(c_dir + "bin/unittests")

if not os.path.isdir(c_dir + "bin/integration"):
    os.mkdir(c_dir + "bin/integration")

#libs = "base math data_structure util vienna secondary_structure eternabot structure motif resources motif_data_structures thermo_fluctuation motif_state_search sequence_optimizer instances"
libs = "base math data_structure util vienna secondary_structure eternabot structure motif motif_tools resources motif_data_structures thermo_fluctuation motif_state_search sequence_optimizer"
lib_paths = libs.split()

depends = {
    'base' : '',
    'math' : 'base',
    'data_structure' : 'base',
    'util' : 'math ${SQLITE3_LIBRARY}',
    'vienna' : 'base',
    'secondary_structure' : 'util',
    'eternabot' : 'vienna secondary_structure',
    'structure' : 'util',
    'motif' : 'structure secondary_structure',
    'motif_tools' : 'motif',
    'resources' : 'motif',
    'motif_data_structures' : 'resources data_structure',
    'thermo_fluctuation' : 'motif_data_structures',
    'motif_state_search' : 'motif_data_structures',
    'sequence_optimizer' : 'motif_data_structures eternabot',
    'all_libs' : 'motif_tools thermo_fluctuation motif_state_search sequence_optimizer',
}

f = open("build.cmake", "w")

for p in lib_paths:

    f.write("set("+p+"_files\n")
    for root, dirnames, filenames in os.walk('../../src/'+p):
        for filename in filenames:
            if filename[-2:] == 'cc' or filename[-3:] == 'cpp':
                f.write("\t"+os.path.join(root, filename)+"\n")
    f.write(")\n")

    f.write("add_library(%s SHARED ${%s_files})\n" % (p, p))
    f.write("target_link_libraries(%s %s)\n\n" % (p, depends[p]))

    unittest_apps = get_unittests_for_dir(p)
    for unit_app in unittest_apps:
        fname = util.filename(unit_app)
        spl = fname.split(".")
        prog_name = spl[0]

        symlink = prog_name+"_symlink"
        f.write("add_executable("+ prog_name + " " + unit_app +")\n")
        f.write("target_link_libraries("+prog_name+" %s)\n" % p)
        f.write("add_custom_target("+symlink+" ALL)\n")
        f.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND ./symlink.py "+prog_name + " " + BIN_PATH + "unittests/"+prog_name+")\n")
        f.write("\n\n")

    integration_apps = get_integration_tests_for_dir(p)
    for inte_app in integration_apps:
        fname = util.filename(inte_app)
        spl = fname.split(".")
        prog_name = spl[0]

        symlink = prog_name+"_symlink"
        f.write("add_executable("+ prog_name + " " + inte_app +")\n")
        f.write("target_link_libraries("+prog_name+" %s)\n" % p)
        f.write("add_custom_target("+symlink+" ALL)\n")
        f.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND ./symlink.py "+prog_name + " " + BIN_PATH +  "integration/"+prog_name+")\n")
        f.write("\n\n")

    #compile unittests for this lib

    #f.write("target_link_libraries ( " + )

f.write("add_library(all_libs SHARED ../../src/main.cpp)\n")
f.write("target_link_libraries(all_libs %s)\n\n" % (depends['all_libs']))

for p in ['apps']:

    unittest_apps = get_unittests_for_dir(p)
    for unit_app in unittest_apps:
        fname = util.filename(unit_app)
        spl = fname.split(".")
        prog_name = spl[0]
        spl2 = prog_name.split("_")
        spl2.pop()
        act_name = "_".join(spl2)
        act_prog_dir = "../../apps/"+act_name+"/"
        #print act_prog_dir
        cc_files = get_cc_files_in_dir(act_prog_dir)
        final_cc_files = []
        for cc_file in cc_files:
            f_cc = open(act_prog_dir+cc_file)
            lines = f_cc.readlines()
            f_cc.close()

            found = 0
            for l in lines:
                if l[0:8] == "int main":
                    found = 1
                    break
            if not found:
                final_cc_files.append(cc_file)

        symlink = prog_name+"_symlink"
        f.write("add_executable("+prog_name + " " + unit_app)
        for fname in cpp_files:
            f.write(act_prog_dir + "/" + fname + " ")


        f.write("target_link_libraries("+prog_name+" %s)\n" % 'all_libs')
        f.write("add_custom_target("+symlink+" ALL)\n")
        f.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND ./symlink.py "+prog_name + " " + BIN_PATH + "unittests/"+prog_name+")\n")
        f.write("\n\n")

f.close()

write_application_cmake_file()