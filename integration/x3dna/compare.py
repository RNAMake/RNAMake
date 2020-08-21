import pandas as pd
import re

class WCPairing:
    def __init__(self,**kwargs):
        self.index1 = 0
        self.index2 = 0
        self.x = 0.
        self.y = 0.
        self.z = 0.

        for variable, value in kwargs.items():
            setattr(self,variable,value)

def make_pairings_wc_only(raw_string):

    if not isinstance(raw_string,str):
        return dict()

    pairs = dict()

    for pair_tk in [tk for tk in raw_string.split('_') if len(tk)] :
        tokens = pair_tk.split('|')
        if len(tokens) < 5: print(tokens);
        #if (wc_only and tokens[0].find("cW-W") != -1) or (not wc_only):
            #print(tokens)
        if tokens[0].find("cW-W") != -1:
            coords = [float(val) for val in tokens[3].split(' ') if len(val) > 0]
            assert len(coords) == 3
            index1 = int(tokens[1])
            index2 = int(tokens[2])

            key = (min(index1,index2),max(index1,index2))

            if key not in pairs:
                pairs[key] = WCPairing(**{
                    "index1" : index1,
                    "index2" : index2,
                    "x" : coords[0],
                    "y" : coords[1],
                    "z" : coords[2],
                    })
    return pairs

def make_all_pairings(raw_string):

    if not isinstance(raw_string,str):
        return dict()

    pairs = dict()

    for pair_tk in [tk for tk in raw_string.split('_') if len(tk)] :
        tokens = pair_tk.split('|')
        if len(tokens) < 5: print(tokens);
        #if (wc_only and tokens[0].find("cW-W") != -1) or (not wc_only):
            #print(tokens)
        if tokens[0].find("cW-W") != -1:
            coords = [float(val) for val in tokens[3].split(' ') if len(val) > 0]
            assert len(coords) == 3
            index1 = int(tokens[1])
            index2 = int(tokens[2])

            key = (min(index1,index2),max(index1,index2))

            if key not in pairs:
                pairs[key] = WCPairing(**{
                    "index1" : index1,
                    "index2" : index2,
                    "x" : coords[0],
                    "y" : coords[1],
                    "z" : coords[2],
                    })
    return pairs

def main():
    df = pd.read_csv("comp2.csv")
    ### things to do. kick out the straight matches.
    json_has_orig = 0
    mismatch = 0
    for i,row in df.iterrows():
        orig_dict = make_pairings_wc_only(row["unresolved_orig_description"])
        json_dict = make_pairings_wc_only(row["unresolved_json_description"])
        valid = re.compile(r"[A-Z0-9]{4}")

        # check for valid PBD code
        if not valid.match(row['pdb_code']):
            raise TypeError(f"Invalid PDB code {row['pdb_code']}")


        # key comparison
        orig_keys = set(orig_dict.keys())
        json_keys = set(json_dict.keys())
        # first, check that all the original keys are in the new keys
        for key in orig_keys:
            if key not in json_keys:
                alt_dict = make_all_pairings(row['unresolved_json_description'])
                alt_keys = set(alt_dict.keys())
                if key not in alt_keys:
                    mismatch += 1
                    break
        else:
            json_has_orig += 1

    print(json_has_orig,mismatch)



if __name__ == "__main__":
    main()
