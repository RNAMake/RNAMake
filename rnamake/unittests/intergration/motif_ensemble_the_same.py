import unittest

from rnamake import resource_manager as rm
from rnamake import sqlite_library, basepair

class MotifEnsembleTests(unittest.TestCase):

    def _generate_profile(self):
        f = open("correct_mes.dat", "w'")

        me_lib = sqlite_library.MotifEnsembleSqliteLibrary("bp_steps")
        bps = "CG,GC,AU,UA,GU,UG".split(",")
        end_ids = []

        for bp1 in bps:
            for bp2 in bps:
                end_ids.append(bp1[0]+bp2[0]+"_LL_"+bp2[1]+bp1[1]+"_RR")

        for end_id in end_ids:
            me = me_lib.get(name=end_id)
            f.write(me.id + "|")
            for mem in me.members:
                f.write(mem.motif.ends[1].state().to_str() + "|")
            f.write("\n")

        #me = me_lib.get(name="GG_LL_CC_RR")
        #print len(me.members)

    def test_motif_state_members(self):
        #self._generate_profile()
        f = open("correct_mes.dat")
        lines = f.readlines()
        f.close()

        count = 0
        me_lib = sqlite_library.MotifEnsembleSqliteLibrary("bp_steps")
        for l in lines:
            spl = l.split("|")
            name = spl.pop(0)
            me = me_lib.get(name=name)
            for i in range(0, len(spl)-1):
                bp_state = basepair.str_to_basepairstate(spl[i])
                mem = me.members[i]
                diff = bp_state.diff(mem.motif.ends[1].state())
                #self.failIf(diff > 0.1)
                if diff > 0.1:
                    count += 1
        print count





def main():
    unittest.main()

if __name__ == '__main__':
    main()
