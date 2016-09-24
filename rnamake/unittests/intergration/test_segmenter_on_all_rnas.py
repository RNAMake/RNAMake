import unittest
import os
import warnings
import rnamake.segmenter
import rnamake.pose
import rnamake.settings
import rnamake.pose_factory as pf
import rnamake.motif_type as motif_type

from rnamake import exceptions, user_warnings



class SegmenterUnittest(unittest.TestCase):

    def test_apply(self):

        # strange ends 1ZX7 354D
        # possibly missing end to bad res: 1ZH5 1ZBI 1Z7F 1YZD 1YZ9
        # looks like circle 1JBS
        # different format: 2GQ6 2NZ4 2VQE 2Z75 3ADD 3FS0
        skip_structs = "354D 1ZX7 1YTU 1YVP 1YZ9 1YZD 1Z7F 1ZBI 1ZH5 1JBS 3SLQ 1B2M 1G2J 1GTF 1GTN 1H2C 1J6S 1KFO 1KNZ 1L3Z 1MDG 1OSU 1P79 1SDS 1UTD 1UTF 1UVI 1UVJ 1UVK 1UVL 1UVM 1ZBL 1ZEV 255D 283D 2A8V 2DB3 2EZ6 2G4B 2GQ6 2J0S 2JEA 2JLU 2JLX 2JLZ 2NZ4 2PO1 2Q1O 2Q1R 2Q66 2VNU 2VOD 2VON 2VOO 2VQE 2X1A 2X1F 2XNR 2XZO 2Z75 3ADD 3AHU 3B0U 3BNT 3BSX 3CZW 3D2S 3D2V 3ER9 3FHT 3FS0 3GIB 3GPQ 3HGA 3HSB 3I5X 3I5Y 3K49 3K5Q 3K5Y 3K5Z 3K61 3K62 3K64 3L26 3M7N 3MJ0 3ND3 3NJ6 3NMR 3NNA 3NNC 3O7V 3O8C".split()
        too_large = "1S72".split()
        path = "/Users/josephyesselman/projects/REDESIGN/resources/non-redundant-rnas"
        dirs = list(os.walk(path))[0][1]
        for d in dirs:
            if d in skip_structs:
                continue
            if d in too_large:
                continue
            #if '3BNO' != d:
            #    continue
            try:
                with warnings.catch_warnings(record=True) as w:
                    p = pf.factory.pose_from_file(path+"/"+d)
            except exceptions.X3dnaException:
                continue

            twoways = p.motifs(motif_type.TWOWAY)
            for i, t in enumerate(twoways):
                if len(t.ends) != 2:
                    continue
                if len(t.chains()) != 2:
                    continue
                #t.to_pdb("test.pdb")
                s = rnamake.segmenter.Segmenter()
                segments = s.apply(p, t.ends)
                self.failUnless(len(t.residues()) == len(segments.removed.residues()))

            """nways = p.motifs(motif_type.NWAY)
            for i, t in enumerate(twoways):
                s = rnamake.segmenter.Segmenter()
                segments = s.apply(p, t.ends)
                self.failUnless(len(t.residues()) == len(segments.removed.residues()))"""









def main():
    unittest.main()

if __name__ == '__main__':
    main()
