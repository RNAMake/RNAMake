import os
import glob
import fnmatch

from rnamake import util

# make sure these directories exist
if not os.path.isdir("../../bin"):
    os.mkdir("../../bin")

if not os.path.isdir("../../bin/unittests"):
    os.mkdir("../../bin/unittests")

libs = "base math data_structure util vienna secondary_structure eternabot structure motif resources motif_data_structures thermo_fluctuation motif_state_search sequence_optimizer instances"
#libs = "base math data_structure util secondary_structure structure"
lib_paths = libs.split()

for p in lib_paths:
    f = open(p+".cmake", "w")
    f.write("set("+p+"_files\n")
    for root, dirnames, filenames in os.walk('../../src/'+p):
        for filename in filenames:
            if filename[-2:] == 'cc' or filename[-3:] == 'cpp':
                f.write("\t"+os.path.join(root, filename)+"\n")
    f.write(")\n")
    f.close()


matches = []
unittest_apps = []
for root, dirnames, filenames in os.walk('../../unittests_new'):
    for filename in fnmatch.filter(filenames, '*.c*'):
        path = os.path.join(root, filename)
        if filename == "all.cpp":
            unittest_apps.append((os.path.join(root, filename)))
            continue

        f = open(path)
        fail = 0
        for l in f.readlines():
            if len(l) < 8:
                continue
            if l[:9] == 'TEST_CASE':
                unittest_apps.append((os.path.join(root, filename)))
                fail = 1
                break
        if fail:
            continue
        matches.append(os.path.join(root, filename))


f = open("unittests.cmake", "w")
f.write("set(unittests_files\n")
for m in matches:
    f.write("\t"+m+"\n")
f.write(")\n")


f = open("libraries.cmake", "w")
f.write("set(libraries ")
f.write(libs + " unittests)\n")
f.close()

fsum = open("apps.cmake", "w")

f = open("apps.txt")
for l in f.readlines():
    if l[0] == "#":
        continue
    spl = l.split()
    symlink = spl[0]+"_symlink"
    fsum.write("add_executable("+l.rstrip()+")\n")
    fsum.write("target_link_libraries("+spl[0]+" ${link_libraries})\n")
    fsum.write("add_custom_target("+symlink+" ALL)\n")
    #fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND cmake -E create_symlink ../../bin/"+spl[0]+" "+ spl[0]+")\n")
    fsum.write("\n\n")

fsum.close()
f.close()

fsum = open("unittest_apps.cmake", "w")

f = open("unittests.txt")

"""for l in f.readlines():
    spl = l.split()
    if len(spl) < 2:
        continue
    symlink = spl[0]+"_symlink"
    fsum.write("add_executable("+l.rstrip()+")\n")
    fsum.write("target_link_libraries("+spl[0]+" ${link_libraries})\n")
    fsum.write("add_custom_target("+symlink+" ALL)\n")
    #fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND  -E create_symlink ../../bin/unittests/"+spl[0]+" "+ spl[0]+")\n")
    fsum.write("\n\n")

fsum.close()"""

for path in unittest_apps:
    fname = util.filename(path)
    if fname == "all_tests.cc" or fname == "main.cpp":
        continue

    spl = fname.split(".")
    if spl[0][-3:] == "app":
        prog_name = spl[0][:-4]
    elif fname == "main.cc":
        prog_name = "test"
    else:
        prog_name = spl[0]

    symlink = prog_name+"_symlink"
    fsum.write("add_executable("+ prog_name + " " + path +")\n")
    fsum.write("target_link_libraries("+prog_name+" ${link_libraries})\n")
    fsum.write("add_custom_target("+symlink+" ALL)\n")
    fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND ./symlink.py "+prog_name + " ../../bin/unittests/"+spl[0]+" "+ spl[0]+")\n")
    fsum.write("\n\n")


f.close()

