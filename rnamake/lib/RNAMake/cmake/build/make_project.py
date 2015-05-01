import os
import glob
import fnmatch

libs = "base math util structure motif"
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
