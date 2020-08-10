//
// Created by Joseph Yesselman on 2019-04-19.
//

#ifndef RNAMAKE_NEW_TECTO_GRAPH_SETUP_H
#define RNAMAKE_NEW_TECTO_GRAPH_SETUP_H

#include <base/log.h>
#include <secondary_structure/pose.h>
#include <motif/motif.h>
#include <motif_data_structure/motif_graph.h>
#include <data_structure/graph_base.h>
#include <resources/resource_manager.h>

//motif::MotifOPs
//get_motifs_from_secondary_structure(
//        secondary_structure::Pose const &);
motif::MotifOPs
get_motifs_from_secondary_structure(
        secondary_structure::Pose const & p) {
    auto & rm = resources::Manager::instance();
    auto motifs = motif::MotifOPs();
    auto motif = motif::MotifOP(nullptr);
    auto start = 0; // keeps track to see if we found the tetraloop-receptor to start recording motifs

    // TODO this does not support 3-way junctions and any complex topology
    for(auto const & m : p.motifs()) {
        if(m->mtype() == util::MotifType::TWOWAY && !start) {
            LOG_VERBOSE << "found tetraloop-receptor -> sequence: " + m->sequence() + " structure: " + m->dot_bracket();
            LOG_VERBOSE << "starting to record motifs ------------------------------------------------";
            start = 1;
            continue;
        }
        if(m->mtype() == util::MotifType::HAIRPIN) {
            LOG_VERBOSE << "found loop -> sequence: " + m->sequence() + " structure: " + m->dot_bracket();
            LOG_VERBOSE << "ending to record motifs   ------------------------------------------------";
            break;
        }
        if(!start) { continue; }
        if(m->mtype() == util::MotifType::HELIX) {
            LOG_VERBOSE << "found base-pair step -> sequence: " + m->sequence();
            try {
                motif = rm.bp_step(m->end_ids()[0]);
            }
            catch(resources::SqliteLibraryException) {
                LOG_ERROR << "could not find base-pair step -> sequence: " + m->sequence() + " is it Watson-Crick?";
                exit(0);
            }
        }
        else if(m->mtype() == util::MotifType::TWOWAY) {
            LOG_VERBOSE << "found twoway junction -> sequence: " + m->sequence() + " structure: " + m->dot_bracket();
            try {
                motif = rm.motif("", m->end_ids()[0]);
            }
            catch(resources::ResourceManagerException) {
                LOG_ERROR << "could not find twoway junction -> sequence: " + m->sequence() + " structure: " + m->dot_bracket();
                exit(0);
            }
        }
        else {
            LOG_ERROR << "motif not supported in tecto simulations -> sequence: " + m->sequence() + " structure: " + m->dot_bracket();
            exit(0);
        }
        motifs.push_back(motif);
    }

    return motifs;
}

struct TectoGraphInfo {
    inline
    TectoGraphInfo(
            motif_data_structure::MotifGraphOP n_mg,
            data_structure::NodeIndexandEdge const & n_start,
            data_structure::NodeIndexandEdge const & n_end):
            mg(n_mg),
            start(n_start),
            end(n_end) {}

    motif_data_structure::MotifGraphOP mg;
    data_structure::NodeIndexandEdge start, end;
};

typedef std::shared_ptr<TectoGraphInfo> TectoGraphInfoOP;

class TectoGraphSetup {
public:
    TectoGraphSetup() {}

    virtual
    ~TectoGraphSetup() {}

    virtual
    TectoGraphSetup *
    clone() const = 0;


public:
    virtual
    TectoGraphInfoOP
    get_graph(
            secondary_structure::Pose const &,
            secondary_structure::Pose const &) = 0;

};

typedef std::shared_ptr<TectoGraphSetup> TectoGraphSetupOP;

class DefaultSetup : public TectoGraphSetup {
public:
    DefaultSetup() : TectoGraphSetup() {}

    ~DefaultSetup() {}

    TectoGraphSetup *
    clone() const override { return new DefaultSetup(*this); };

public:
    TectoGraphInfoOP
    get_graph(
            secondary_structure::Pose const & flow_p,
            secondary_structure::Pose const & chip_p) override {
        LOG_VERBOSE << "processing sequence of flow RNA";
        auto flow_motifs = get_motifs_from_secondary_structure(flow_p);
        LOG_VERBOSE << "processing sequence of chip RNA";
        auto chip_motifs = get_motifs_from_secondary_structure(chip_p);

        auto & rm = resources::Manager::instance();
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        mg->set_option_value("sterics", false);
        mg->add_motif(rm.bp_step("GG_LL_CC_RR"));
        mg->add_motif(rm.motif("GGAA_tetraloop", "", "A14-A15"));
        mg->add_motif(flow_motifs[1], 1,  "A7-A22");
        //std::cout << mg->get_node(1)->data()->ends()[1]->name() << std::endl;
        for(int i = 2; i < flow_motifs.size(); i++) { mg->add_motif(flow_motifs[i]); }
        mg->add_motif(rm.motif("GAAA_tetraloop", "", "A149-A154"));
        mg->add_motif(chip_motifs[1], -1, "A222-A251");
        for(int i = 2; i < chip_motifs.size(); i++) { mg->add_motif(chip_motifs[i]); }
        //mg->add_connection(1, mg->last_node()->index(), "A1-A6", mg->last_node()->data()->end_name(1));
        auto start = data_structure::NodeIndexandEdge{1, mg->get_node(1)->data()->get_end_index("A1-A6") };
        auto end   = data_structure::NodeIndexandEdge{mg->last_node()->index(), 1};

        return std::make_shared<TectoGraphInfo>(mg, start, end);

    }

};


#endif //RNAMAKE_NEW_TECTO_GRAPH_SETUP_H
