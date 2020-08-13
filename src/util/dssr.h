#ifndef __DSSR_H__
#define __DSSR_H__

#include <vector>
#include <limits>

#include <math/quaternion.h>
#include <math/xyz_vector.h>
#include <base/types.h>
#include <base/sys_interface.h>

namespace util {


double
get_double(const nlohmann::json& , const String&);

int
get_int(const nlohmann::json& , const String& );

char
get_char(const nlohmann::json& , const String&);

String
get_string(const nlohmann::json& , const String&);

math::Point
get_point(const nlohmann::json& , const String&);

math::Quaternion
get_quaternion(const nlohmann::json&, const String&);

math::Matrix
get_matrix(const nlohmann::json&);

struct DssrNt {
    math::Point C5prime_xyz{-1.,-1.,-1.};
    double Dp{-1.};
    math::Point P_xyz{-1.,-1.,-1.};
    double alpha{-1.};
    double amplitude{-1.};
    String bb_type{"NA"}; // TODO make into an enum
    double beta{-1.};
    String bin	{"NA"};
    char chain_name{' '};
    double chi{-1.};
    String cluster{"NA"};
    char dbn{' '};
    double delta{-1.};
    double epsilon{-1.};
    double epsilon_zeta{-1.};
    double eta{-1.};
    double eta_base{-1.};
    double eta_prime{-1.};
    double filter_rmsd{-1.};
    char form{' '};
    math::Point frame_origin{-1.,-1.,-1.};
    math::Quaternion frame_quaternion{-1.,-1.,-1.,-1.};
    double frame_rmsd {-1.}; 
    math::Matrix ref_frame {-1.,-1.,-1., -1,-1.,-1., -1.,-1.,-1.};
    double gamma {-1.};
    String glyco_bond {"NA"};
    int index{-1};
    int index_chain{-1};
    char nt_code{' '};
    String nt_id{"NA"};
    char nt_name{' '};
    int nt_resnum{-1};
    double phase_angle{-1.};
    String  puckering{"NA"};
    double splay_angle{-1.};
    double splay_distance{-1.};
    double splay_ratio{-1.};
    double ssZp{-1.};
    String sugar_class{"NA"};
    double suiteness{-1.};
    String summary{"NA"};
    double theta{-1.};
    double theta_null{-1.};
    double theta_prime{-1.};
    double v0{-1.};
    double v1{-1.};
    double v2{-1.};
    double v3{-1.};
    double v4{-1.};
    double zeta{-1.};
};

struct DssrPair {


};

using DssrNts = std::vector<DssrNt>;
using DssrPairs = std::vector<DssrPair>;

DssrNts
get_nts(const nlohmann::json& );

DssrPairs
get_pairs(const nlohmann::json&);

void
get_elements(
        String const & ,
        DssrNts & ,
        DssrPairs &

        ) ;
    //from typing import List, Dict
//
//
//class DSSR_NT (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.index_chain : int = None
//        self.chain_name : str = None
//        self.nt_resnum : int = None
//        self.nt_name : str = None
//        self.nt_code : str = None
//        self.nt_id : str = None
//        self.dbn : str = None
//        self.summary : str = None
//        self.alpha = None
//        self.beta = None
//        self.gamma : float = None
//        self.delta : float = None
//        self.epsilon : float = None
//        self.zeta : float = None
//        self.epsilon_zeta : float = None
//        self.bb_type : str = None
//        self.chi : float = None
//        self.glyco_bond : str = None
//        self.C5prime_xyz : list = None
//        self.P_xyz : list = None
//        self.form : str = None
//        self.ssZp : float = None
//        self.Dp : float = None
//        self.splay_angle : float = None
//        self.splay_distance : float = None
//        self.splay_ratio : float = None
//        self.eta = None
//        self.theta = None
//        self.eta_prime = None
//        self.theta_prime = None
//        self.eta_base = None
//        self.theta_base = None
//        self.v0 : float = None
//        self.v1 : float = None
//        self.v2 : float = None
//        self.v3 : float = None
//        self.v4 : float = None
//        self.amplitude : float = None
//        self.phase_angle : float = None
//        self.puckering : str = None
//        self.sugar_class : str = None
//        self.bin : str = None
//        self.cluster : str = None
//        self.suiteness : float = None
//        self.filter_rmsd : float = None
//        self.frame : dict = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_PAIR (object):
//    def __init__(self, nt1 : DSSR_NT, nt2 : DSSR_NT, **kwargs):
//        self.nt1 = nt1
//        self.nt2 = nt2
//        self.index : int = None
//        self.bp : str = None
//        self.name : str = None
//        self.Saenger : str = None
//        self.LW : str = None
//        self.DSSR : str = None
//        self.chi1 : float = None
//        self.conf1 : str = None
//        self.pucker1 : str = None
//        self.lambda1 : float = None
//        self.chi2 : float = None
//        self.conf2 : str = None
//        self.pucker2 : str = None
//        self.lambda2 : float = None
//        self.C1C1_dist : float = None
//        self.N1N9_dist : float = None
//        self.C6C8_dist : float = None
//        self.CNNC_torsion : float = None
//        self.hbonds_num : int = None
//        self.hbonds_desc : str = None
//        self.interBase_angle : float = None
//        self.planarity : float = None
//        self.simple_Shear : float = None
//        self.simple_Stretch : float = None
//        self.simple_Buckle : float = None
//        self.simple_Propeller : float = None
//        self.bp_params : list = None
//        self.frame : dict = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_HAIRPIN (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.type : str = None
//        self.bridging_nts : list = None
//        self.stem_indices : list = None
//        self.summary : str = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//        self.num_stems : int = None
//        self.bridges : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_HELIX (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.num_stems : int = None
//        self.strand1 : str = None
//        self.strand2 : str = None
//        self.bp_type : str = None
//        self.helix_form : str = None
//        self.helical_rise : float = None
//        self.helical_rise_std : float = None
//        self.helical_radius : float = None
//        self.helical_radius_std : float = None
//        self.helical_axis : list = None
//        self.point1 : list = None
//        self.point2 : list = None
//        self.num_pairs : int = None
//        self.pairs : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_STEM (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.helix_index : int = None
//        self.strand1 : str = None
//        self.strand2 : str = None
//        self.bp_type : str = None
//        self.helix_form : str = None
//        self.helical_rise : float = None
//        self.helical_rise_std : float = None
//        self.helical_radius : float = None
//        self.helical_radius_std : float = None
//        self.helical_axis : list = None
//        self.point1 : list = None
//        self.point2 : list = None
//        self.num_pairs : int = None
//        self.pairs : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_ILOOP (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.type : str = None
//        self.bridging_nts : list = None
//        self.stem_indices : list = None
//        self.summary : str = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//        self.num_stems : int = None
//        self.bridges : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_JUNCTION (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.type : str = None
//        self.bridging_nts : list = None
//        self.stem_indices : list = None
//        self.summary : str = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//        self.num_stems : int = None
//        self.bridges : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_SINGLE_STRAND (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_KISSING_LOOP (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.stem_index : int = None
//        self.hairpin_indices : list = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_AMINOR (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.type : str = None
//        self.desc_short : str = None
//        self.desc_long : str = None
//        self.A_nt1 : dict = None
//        self.A_nt2 : dict = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_RIBOSE_ZIPPER (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_PSEUDOKNOT (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.desc : str = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_HBOND (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.atom1_serNum : int = None
//        self.atom2_serNum : int = None
//        self.donAcc_type : str = None
//        self.distance : float = None
//        self.atom1_id : str = None
//        self.atom2_id : str = None
//        self.atom_pair : str = None
//        self.residue_pair : str = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
//class DSSR_SPLAY_UNITS (object):
//    def __init__(self, **kwargs):
//        self.index : int = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//
//        for key, value in kwargs.items():
//            setattr(self, key, value)
//
//
}
#endif// __DSSR_H__
