import exceptions


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
|   SSTRAND         | 10    | Single stranded RNA                         |
+-------------------+-------+---------------------------------------------+
|   TCONTACT        | 11    | General tertiary contact                    |
+-------------------+-------+---------------------------------------------+
|   UNKNOWN         | 99    | Motif with unassigned motif type            |
+-------------------+-------+---------------------------------------------+
|   ALL             | 999   | Place holder type for iteration of motifs   |
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
    TWOWAY            : 'TWOWAY',
    NWAY              : 'NWAY',
    HAIRPIN           : 'HAIRPIN',
    TCONTACT          : 'TCONTACT',
    TCONTACT_HP_HP    : 'TCONTACT_HP_HP',
    TCONTACT_H_HP     : 'TCONTACT_H_HP',
    TCONTACT_H_H      : 'TCONTACT_H_H',
    TWOWAY_SEGMENTS   : 'TWOWAY_SEGMENTS',
    HELIX             : 'HELIX',
    UNKNOWN           : 'UNKNOWN',
    ALL               : 'ALL'
}

str_to_type_dict = {
   'TWOWAY'           : TWOWAY,
   'NWAY'             : NWAY,
   'HAIRPIN'          : HAIRPIN,
   'TCONTACT_HP_HP'   : TCONTACT_HP_HP,
   'TCONTACT_H_HP'    : TCONTACT_H_HP,
   'TCONTACT_H_H'     : TCONTACT_H_H,
   'TWOWAY_SEGMENTS'  : TWOWAY_SEGMENTS,
   'HELIX'            : HELIX,
   'UNKNOWN'          : UNKNOWN,
   'ALL'              : ALL
}


def type_to_str(mtype):
    """
    convery motif_type enum value into a string for each printing

    :param mtype: motif type value
    :type mtype: motif_type
    :return: return string name of motif_type from type_to_str_dict
    :rtype: string

    :examples:

    ..  code-block:: python

        >>> import motif_type
        >>> motif_type.type_to_str(motif_type.TWOWAY)
        TWOWAY

    """

    if mtype in type_to_str_dict:
        return type_to_str_dict[mtype]
    else:
        raise exceptions.MotifTypeException("MotifType not recognized: " + mtype)


def str_to_type(type_name):
    """
    converts string to its corresponding motif_type

    :param type_name: name of motif_type enum
    :return: corresponding motif_type enum
    :rtype: motif_type

    ..  code-block:: python

        >>> import motif_type
        >>> motif_type.str_to_type('TWOWAY')
        0

    """

    if type_name in str_to_type_dict:
        return str_to_type_dict[type_name.upper()]
    else:
        raise exceptions.MotifTypeException("motif_type name not recognized")




