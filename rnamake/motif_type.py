


"""
Basic type information for motifs for looking up in databases and libraries

+-------------------+-------+---------------------------------------------+
|  Motif Type       | Value | Description                                 |
+===================+=======+=============================================+
|   TWOWAY          | 0     | Twoway junctions                            |
+-------------------+-------+---------------------------------------------+
|   NWAY            | 1     | Junctions with more then two ends           |
+-------------------+-------+---------------------------------------------+
|   HAIRPIN         | 2     | Motif with only one end                     |
+-------------------+-------+---------------------------------------------+
|   TCONTACT_HP_HP  | 3     | Tertiary contact between two hairpins,      |
|                   |       | two ends                                    |
+-------------------+-------+---------------------------------------------+
|   TCONTACT_H_HP   | 4     | Tertiary contact between one hairpin and    |
|                   |       | one helix, three ends                       |
+-------------------+-------+---------------------------------------------+
|   TCONTACT_H_H    | 5     | Tertiary contact between two helices,       |
|                   |       | four ends                                   |
+-------------------+-------+---------------------------------------------+
|   TWOWAY_SEGMENTS | 8     | Segments of multiple motifs with two ends   |
+-------------------+-------+---------------------------------------------+
|   HELIX           | 9     | RNA with only GC/AU pairs                   |
+-------------------+-------+---------------------------------------------+
"""

TWOWAY            = 0
NWAY              = 1
HAIRPIN           = 2
TCONTACT_HP_HP    = 3
TCONTACT_H_HP     = 4
TCONTACT_H_H      = 5
TWOWAY_SEGMENTS   = 8
HELIX             = 9
SSTRAND           = 10
TCONTACT          = 11
UNKNOWN           = 99
ALL               = 999

type_to_str_dict = {
    TWOWAY           : 'TWOWAY',
    NWAY             : 'NWAY',
    HAIRPIN          : 'HAIRPIN',
    TCONTACT         : 'TCONTACT',
    TCONTACT_HP_HP   : 'TCONTACT_HP_HP',
    TCONTACT_H_HP    : 'TCONTACT_H_HP',
    TCONTACT_H_H     : 'TCONTACT_H_H',
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
        raise ValueError("motif_type name not recognized")

def is_valid_motiftype(mtype):
    type_to_str(mtype)



