import os
import glob
import fnmatch

libs = "base math util structure motif motif_tree_state resources motif_assembly"
lib_paths = libs.split()

for p in lib_paths:
    files = glob.glob("../../src/"+p+"/*.cc")
    f = open(p+".cmake", "w")
    f.write("set("+p+"_files\n")
    for file in files:
        f.write("\t"+file+"\n")
    f.write(")\n")

matches = []
progs = "main.cc all_tests.cc".split()
for root, dirnames, filenames in os.walk('../../unittests'):
    for filename in fnmatch.filter(filenames, '*.cc'):
        if filename in progs:
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
    spl = l.split()
    symlink = spl[0]+"_symlink"
    fsum.write("add_executable("+l.rstrip()+")\n")
    fsum.write("target_link_libraries("+spl[0]+" ${link_libraries})\n")
    fsum.write("add_custom_target("+symlink+" ALL)\n")
    fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND cmake -E create_symlink ../../bin/"+spl[0]+" "+ spl[0]+")\n")
    fsum.write("\n\n")

fsum.close()
f.close()

fsum = open("unittest_apps.cmake", "w")

f = open("unittests.txt")
for l in f.readlines():
    spl = l.split()
    symlink = spl[0]+"_symlink"
    fsum.write("add_executable("+l.rstrip()+")\n")
    fsum.write("target_link_libraries("+spl[0]+" ${link_libraries})\n")
    fsum.write("add_custom_target("+symlink+" ALL)\n")
    fsum.write("add_custom_command(TARGET "+symlink+" POST_BUILD COMMAND cmake -E create_symlink ../../bin/unittests/"+spl[0]+" "+ spl[0]+")\n")
    fsum.write("\n\n")

fsum.close()
f.close()

