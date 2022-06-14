#ifndef __DSSR_H__
#define __DSSR_H__

#include <vector>
#include <limits>
#include <iostream>
#include <limits>

#include <math/quaternion.h>
#include <math/vector_3.hpp>
#include <base/types.hpp>
#include <math/matrix_3x3.hpp>
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

math::Vector3
get_point(const nlohmann::json& , const String&);

math::Quaternion
get_quaternion(const nlohmann::json&, const String&);

math::Matrix3x3
get_matrix(const nlohmann::json&);

Reals
get_reals(const nlohmann::json&, const String&);

Indexes
get_ints(const nlohmann::json&, const String&);


struct DssrElem {

    virtual
    bool
    valid() const  = 0;
};

struct DssrNt : public DssrElem {
    math::Vector3 C5prime_xyz{-1.,-1.,-1.};
    double Dp{-1.};
    math::Vector3 P_xyz{-1.,-1.,-1.};
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
    math::Vector3 frame_origin{-1.,-1.,-1.};
    math::Quaternion frame_quaternion{-1.,-1.,-1.,-1.};
    double frame_rmsd {-1.}; 
    math::Matrix3x3 ref_frame {-1.,-1.,-1., -1,-1.,-1., -1.,-1.,-1.};
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
    math::Vector3 frame_origin{-1.,-1.,-1.};
    math::Quaternion frame_quaternion{-1.,-1.,-1.,-1.};
    double frame_rmsd {-1.}; 
    math::Matrix3x3 ref_frame {-1.,-1.,-1., -1,-1.,-1., -1.,-1.,-1.};
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
    Indexes bridging_nts{-1};
    Indexes stem_indices{-1};
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
    math::Vector3 helical_axis{-1.,-1.,-1.};
    math::Vector3 point1{-1.,-1.,-1.};
    math::Vector3 point2{-1.,-1.,-1.};
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
    math::Vector3 helical_axis{-1.,-1.,-1.};
    math::Vector3 point1{-1.,-1.,-1.};
    math::Vector3 point2{-1.,-1.,-1.};
    int num_pairs{-1};
    std::vector<DssrPair> pairs;

    bool 
    valid() const override {
        return false;
    }
};

struct DssrILoop : public DssrElem {
    
    int index{-1};
    String type{"NA"}; 
    Indexes bridging_nts{};
    Indexes stem_indices{};
    String summary{"NA"};
    int num_nts{-1};
    String nts_short{"NA"};
    String nts_long{"NA"};
    int num_stems{-1};

    //        self.bridges : list = None

    bool 
    valid() const override {
        return false;
    }
};

struct DssrJunc : public DssrElem {
    int index{-1};
    String type{"NA"};
    Indexes bridging_nts{};
    Indexes stem_indices{};
    String summary{"NA"};
    int num_nts{-1};
    String nts_short{"NA"};
    String nts_long{"NA"};
    int num_stems{-1};

    //        self.bridges : list = None

    bool 
    valid() const override {
        return false;
    }

};

struct DssrSingStrand : public DssrElem {


    int index{-1};
    int num_nts{-1};
    String nts_short{"NA"};
    String nts_long{"NA"};

    bool 
    valid() const override {
        return false;
    }
};

struct DssrKissingLoop : public DssrElem {

    int index{-1};
    int stem_index{-1};
    Indexes hairpin_indices{};

    bool 
    valid() const override {
        return false;
    }
};

struct DssrAMinor : public DssrElem {
    int index{-1};
    String type{"NA"};
    String desc_short{"NA"};
    String desc_long{"NA"};
//        self.A_nt1 : dict = None TODO these objects... idk what they are yet
//        self.A_nt2 : dict = None
    bool 
    valid() const override {
        return false;
    }
};

struct DssrRiboseZip : public DssrElem {
    int index{-1};
    int num_nts{-1};
    String nts_short{"NA"};
    String nts_long{"NA"};

    bool 
    valid() const override {
        return false;
    }

};

struct DssrPseudoKnot : public DssrElem {

    int index{-1};
    String desc{"NA"};

    bool 
    valid() const override {
        return false;
    }
};

struct DssrHBond : public DssrElem {
    int index{-1};
    int atom1_serNum{-1};
    int atom2_serNum{-1};
    String donAcc_type{"NA"};
    double distance{-1};
    String atom1_id{"NA"};
    String atom2_id{"NA"};
    String atom_pair{"NA"};
    String residue_pair{"NA"};
    
    bool 
    valid() const override {
        return false;
    }

};

struct DssrSplayUnits : public DssrElem {
    int index{-1};
    int num_nts{-1};
    String nts_short{"NA"};
    String nts_long{"NA"};

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
using DssrJuncs = std::vector<DssrJunc>;
using DssrSingStrands = std::vector<DssrSingStrand>;
using DssrKissingLoops = std::vector<DssrKissingLoop>;

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

DssrJuncs
get_juncs(const nlohmann::json&);

void
get_elements(
        String const & ,
        DssrNts & ,
        DssrPairs &,
        DssrHairpins &,
        DssrHelices &,
        DssrStems &,
        DssrILoops &
//        DrrsJuncs &
        ) ;

}
#endif// __DSSR_H__
