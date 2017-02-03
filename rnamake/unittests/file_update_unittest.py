import unittest
from rnamake import settings, file_updater

class FileUpdaterUnittest(unittest.TestCase):
    def test_update_motif_str(self):
        fpath = settings.UNITTEST_PATH+"resources/motifs/start_tar.motif"
        f = open(fpath)
        lines = f.readlines()
        f.close()

        new_s = file_updater.update_motif_str(lines[0])




def main():
    unittest.main()

if __name__ == '__main__':
    main()





