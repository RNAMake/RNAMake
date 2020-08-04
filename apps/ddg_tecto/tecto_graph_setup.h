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

motif::MotifOPs
get_motifs_from_secondary_structure(
        secondary_structure::Pose const &);

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
    clone() const { return new DefaultSetup(*this); };

public:
    TectoGraphInfoOP
    get_graph(
            secondary_structure::Pose const & flow_p,
            secondary_structure::Pose const & chip_p) {
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
