//
// Created by Joseph Yesselman on 2/16/19.
//

#include <util/steric_lookup.hpp>
#include "base/backtrace.hpp"
#include "util/basic_io.hpp"
#include "structure/pdb_parser.h"
#include "resources/resource_manager.h"
#include "dock_motif/dock_motif.h"

DockMotifApp::DockMotifApp() : Application() {}

void
DockMotifApp::setup_options() {

    add_option("scaffold", String(""), OptionType::STRING, true);
    add_option("motif", String(""), OptionType::STRING, true);
    add_option("origin", String(""), OptionType::STRING, false);
    add_option("ligand", String(""), OptionType::STRING, false);
}

void
DockMotifApp::parse_command_line(
        int argc,
        const char ** argv) {

    Application::parse_command_line(argc, argv);

}

void
DockMotifApp::run() {

    auto mf = MotifFactory();
    auto added_helix = mf.added_helix();
    auto motif_path = get_string_option("motif");

    //auto scaffold = RM::instance().motif("scaffold");
    auto scaffold = RM::instance().get_structure(get_string_option("scaffold"), "scaffold");

    RM::instance().add_motif(get_string_option("motif"), "motif", MotifType::HAIRPIN);
    auto motif = RM::instance().motif("motif");
    motif->ends()[0]->flip();
    //auto aligned_m = get_aligned_motif(motif->ends()[0], added_helix->ends()[0], added_helix);

    //auto m_state = motif->get_state();

    auto lookup = StericLookup(1.0, 5.0, 6);
    auto bead_centers = Points();
    for(auto const & b: scaffold->get_beads()) {
        bead_centers.push_back(b.center());
    }
    lookup.add_points(bead_centers);

    //RM::instance().add_motif()

    auto center = Point();
    if(get_string_option("ligand") != "") {
        auto ligand = _parse_ligand_for_center_coords();
        center = ligand->center();
    }

    auto box_size = 12.0;
    auto grid_size = 1.0;
    auto valid_start_points = Points();
    auto current = Point();
    auto clash = 0;
    for(auto i = -box_size; i < box_size; i += grid_size) {
        for(auto j = -box_size; j < box_size; j += grid_size) {
            for(auto k = -box_size; k < box_size; k += grid_size) {
                current = center + Point(i, j, k);
                clash = lookup.clash(current);
                if(clash) { continue; }
                valid_start_points.push_back(current);
            }
        }
    }

    auto last_center = _calc_motif_center(motif);
    for(auto const & p : valid_start_points) {
        motif->move(p - last_center);
        std::cout << last_center << std::endl;
        std::cout << p << std::endl;
        std::cout << p - last_center << std::endl;
        std::cout << _calc_motif_center(motif) << std::endl;
        exit(0);
    }

    //points_to_pdb("grid.pdb", valid_start_points);
}

Point
DockMotifApp::_calc_motif_center(
        MotifOP m ) {
    auto m_center = Point();
    for(auto b : m->beads()) {
        m_center += b.center();
    }
    m_center /= m->beads().size();
    return m_center;
}

ResidueOP
DockMotifApp::_parse_ligand_for_center_coords() {
    auto ligand_pdb_path = get_string_option("ligand");
    auto pdb_parser = PDBParser();
    auto res = pdb_parser.parse(ligand_pdb_path, 0, 0, 1);
    if(res.size() == 0) {
        throw DockMotifAppException(
                ligand_pdb_path + " contains no ligands, make sure you supply one ligand in -ligand path");
    }

    if(res.size() > 1) {
        throw DockMotifAppException(
                ligand_pdb_path + " contains more than one residue, make sure there is only one!");
    }

    return res[0];
}


int
main(
        int argc,
        const char ** argv) {

    //must add this for all apps!
    std::set_terminate(print_backtrace);

    //load extra motifs being used
    //String base_path = base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";

    auto app = DockMotifApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}