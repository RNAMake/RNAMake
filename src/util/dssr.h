#ifndef __DSSR_H__
#define __DSSR_H__

#include <vector>
#include <limits>
#include <iostream>

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

Reals
get_reals(const nlohmann::json&, const String&,int);

Ints
get_ints(const nlohmann::json&, const String&,int);


struct DssrElem {

    virtual
    bool
    valid() const  = 0;
};

struct DssrNt : public DssrElem {
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
    
   
    bool 
    valid() const override {
        //TODO add in the things this class needs to be used by the rest of the programs
        return false;
    }
};

struct DssrPair : public DssrElem {
    double C1C1_dist{-1.};
    double C6C8_dist{-1.};
    double CNNC_torsion{-1};
    String DSSR{"NA"}; //TODO make this an enum	"cW-W"
    String LW{"NA"}; //TODO make this an enum"cWW"
    double N1N9_dist{-1.};//9.305
    String Saenger{"NA"};//	"19-XIX"h
    String bp{"NA"}; //TODO make this an enum	"G-C"
    Reals bp_params{-1.,-1.,-1.,-1.,-1.,-1.}; // TODO make this an array? 
    double chi1 {-1.};	
    double chi2	{-1.}; 
    String conf1{"NA"}; //TODO enum	"anti"
    String conf2{"NA"};	//TODO enum
    math::Point frame_origin{-1.,-1.,-1.};
    math::Quaternion frame_quaternion{-1.,-1.,-1.,-1.};
    double frame_rmsd {-1.}; 
    math::Matrix ref_frame {-1.,-1.,-1., -1,-1.,-1., -1.,-1.,-1.};
    String hbonds_desc{"NA"}; //	"N1(imino)-N3[3.14],N2(amino)-O2(carbonyl)[3.27],O6(carbonyl)-N4(amino)[2.96]"
    int hbonds_num{-1};
    int index{-1};
    double interBase_angle{-1.};
    double lambda1{-1};
    double lambda2{-1};
    String name{"NA"}; // TODO turn this into an enum
    String nt1{"NA"};
    String nt2{"NA"};
    double planarity{-1.};
    String pucker1{"NA"};
    String pucker2{"NA"};
    double simple_Buckle{-1.};
    double simple_Propeller{-1.};
    double simple_Shear{-1.};
    double simple_Stretch{-1.};


    bool 
    valid() const override {
        return false;
    }
};

struct DssrHairpin : public DssrElem {
    //TODO bridges attribute
    int index{-1};
    String nts_long{"NA"};
    String nts_short{"NA"};
    String type{"NA"};
    Ints bridging_nts{-1}; 
    Ints stem_indices{-1}; 
    String summary{"NA"};
    int num_nts{-1};
    int num_stems{-1};

    bool 
    valid() const override {
        return false;
    }
};


struct DssrHelix : public DssrElem {

    int index{-1};
    int num_stems{-1};
    String strand1{"NA"};    
    String strand2{"NA"};    
    String bp_type{"NA"};
    String helix_form{"NA"};
    double helical_rise{-1};    
    double helical_rise_std{-1};    
    double helical_radius{-1};    
    double helical_radius_std{-1};    
    math::Point helical_axis{-1.,-1.,-1.};
    math::Point point1{-1.,-1.,-1.};
    math::Point point2{-1.,-1.,-1.};
    int num_pairs{-1};
    std::vector<DssrPair> pairs;

    bool 
    valid() const override {
        return false;
    }
};

struct DssrStem : public DssrElem {
    int index{-1};
    int helix_index{-1}; 
    String strand1{"NA"};    
    String strand2{"NA"};    
    String bp_type{"NA"};
    String helix_form{"NA"};
    double helical_rise{-1};    
    double helical_rise_std{-1};    
    double helical_radius{-1};    
    double helical_radius_std{-1};    
    math::Point helical_axis{-1.,-1.,-1.};
    math::Point point1{-1.,-1.,-1.};
    math::Point point2{-1.,-1.,-1.};
    int num_pairs{-1};
    std::vector<DssrPair> pairs;

    bool 
    valid() const override {
        return false;
    }
};

struct DssrILoop : public DssrElem {
    int index{-1};
    //        self.type : str = None
//        self.bridging_nts : list = None
//        self.stem_indices : list = None
//        self.summary : str = None
//        self.num_nts : int = None
//        self.nts_short : str = None
//        self.nts_long : str = None
//        self.num_stems : int = None
//        self.bridges : list = None

    bool 
    valid() const override {
        return false;
    }
};

using DssrNts = std::vector<DssrNt>;
using DssrPairs = std::vector<DssrPair>;
using DssrHairpins =  std::vector<DssrHairpin>;
using DssrHelices = std::vector<DssrHelix>;
using DssrStems = std::vector<DssrStem>;
using DssrILoops = std::vector<DssrILoop>;

DssrNts
get_nts(const nlohmann::json& );

DssrPairs
get_pairs(const nlohmann::json&);

DssrHairpins
get_hairpins(const nlohmann::json&);

DssrHelices
get_helices(const nlohmann::json&);

DssrStems
get_stems(const nlohmann::json&);

DssrILoops
get_iloops(const nlohmann::json&);

void
get_elements(
        String const & ,
        DssrNts & ,
        DssrPairs &,
        DssrHairpins &,
        DssrHelices &,
        DssrStems &,
        DssrILoops &
        ) ;

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
