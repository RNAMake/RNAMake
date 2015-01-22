


"""
Basic type information for motifs for looking up in databases and libraries

Attributes
----------
`TWOWAY` : Twoway junctions, an group of interacting residues that are not in watson and crick basepairs and that have only two basepair ends


"""

TWOWAY            = 0
NWAY              = 1
HAIRPIN           = 2
TCONTACT_HP_HP    = 3
TCONTACT_H_HP     = 4
TCONTACT_H_H      = 5
T_T               = 6
T_T_T             = 7
TWOWAY_SEGMENTS   = 8
HELIX             = 9
SSTRAND           = 10
UNKNOWN           = 99
ALL               = 999

type_to_str_dict = {
    TWOWAY           : 'TWOWAY',
    NWAY             : 'NWAY',
    HAIRPIN          : 'HAIRPIN',
    TCONTACT_HP_HP   : 'TCONTACT_HP_HP',
    TCONTACT_H_HP    : 'TCONTACT_H_HP',
    TCONTACT_H_H     : 'TCONTACT_H_H',
    T_T              : '2X_TWOWAY',
    T_T_T            : '3X_TWOWAY',
    TWOWAY_SEGMENTS  : 'TWOWAY_SEGMENTS',
    HELIX            : 'HELIX',
    UNKNOWN          : 'UNKNOWN',
    ALL              : 'ALL'
}

str_to_type_dict = {
   'TWOWAY'                    : TWOWAY,
   'NWAY'                      : NWAY,
   'HAIRPIN'                   : HAIRPIN,
   'TCONTACT_HP_HP'            : TCONTACT_HP_HP,
   'TCONTACT_H_HP'             : TCONTACT_H_HP,
   'TCONTACT_H_H'              : TCONTACT_H_H,
   '2X_TWOWAY'                 : T_T,
   '3X_TWOWAY'                 : T_T_T,
   'TWOWAY_SEGMENTS'           : TWOWAY_SEGMENTS,
   'HELIX'                     : HELIX,
   'UNKNOWN'                   : UNKNOWN,
   'ALL'                       : ALL
}


def type_to_str(mtype):
    if mtype in type_to_str_dict:
        return type_to_str_dict[mtype]
    else:
        raise ValueError("MotifType not recognized: " + mtype)

def str_to_type(type_name):
    if type_name in str_to_type_dict:
        return str_to_type_dict[type_name.upper()]
    else:
        raise ValueError("mtdb_type name not recognized")

def is_valid_motiftype(mtype):
    type_to_str(mtype)



