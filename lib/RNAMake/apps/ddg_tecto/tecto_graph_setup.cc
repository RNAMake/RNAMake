//
// Created by Joseph Yesselman on 2019-04-19.
//

#include <ddg_tecto/tecto_graph_setup.h>
#include <base/log.h>

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

