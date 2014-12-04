import sys

def to_camel_case(name):
    spl = name.split("_")
    name_cc = ""
    for word in spl:
        name_cc +=  word[0].upper() + word[1:]
    return name_cc

def make_unittest(name):
    name_cc = to_camel_case(name)
    f = open(name+"_unittest.py","w")
    string =  "import unittest\nimport rnamake."+name+"\n\n"
    string += "class "+name_cc+"(unittest.TestCase):\n\n"
    string += "    def test_creation(self):\n        pass\n\n"
    string += "def main():\n    unittest.main()\n\n"
    string += "if __name__ == '__main__':\n    main()"
    f.write(string)
    f.close()

def main():
    name = sys.argv[1]
    make_unittest(name)

if __name__ == '__main__':
    main()
