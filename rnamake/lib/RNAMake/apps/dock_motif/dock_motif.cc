//
// Created by Joseph Yesselman on 2/16/19.
//

#include <util/steric_lookup.hpp>
#include <util/monte_carlo.h>
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

    auto scaffold = RM::instance().get_structure(get_string_option("scaffold"), "scaffold");

    RM::instance().add_motif(get_string_option("motif"), "motif", MotifType::HAIRPIN);
    auto motif = RM::instance().motif("motif");
    motif->ends()[0]->flip();

    auto start_ms = motif->get_state();
    _precompute_rotations(start_ms);

    auto screening_lookup = StericLookup(1.0, 5.0, 6);
    lookup_ = StericLookup(1.0, 4.0, 6);
    auto bead_centers = Points();
    for(auto const & b: scaffold->get_beads()) {
        bead_centers.push_back(b.center());
    }
    screening_lookup.add_points(bead_centers);
    lookup_.add_points(bead_centers);

    auto center = Point();
    if(get_string_option("ligand") != "") {
        auto ligand = _parse_ligand_for_center_coords();
        center_ = ligand->center();
    }

    auto ms = _get_starting_state(screening_lookup, center_);
    _search(center_);
    //RM::instance().add_motif()

    auto box_size = 15.0;
    auto grid_size = 1.0;
    auto valid_start_points = Points();
    auto current = Point();
    auto clash = 0;
    for(auto i = -box_size; i < box_size; i += grid_size) {
        for(auto j = -box_size; j < box_size; j += grid_size) {
            for(auto k = -box_size; k < box_size; k += grid_size) {
                current = center + Point(i, j, k);
                clash = screening_lookup.clash(current);
                if(clash) { continue; }
                valid_start_points.push_back(current);
            }
        }
    }



    //points_to_pdb("grid.pdb", valid_start_points);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MotifStateOP
DockMotifApp::_get_starting_state(
        StericLookup & screening_lookup,
        Point const & center) {

    auto ms = rotations_[0];
    auto rng = RandomNumberGenerator();
    auto p = Point();
    int bound = 2;

    ms->move(center - ms->end_states()[0]->d());

    while(1) {
        p = get_random_point(rng, bound);
        ms->move(p);
        if(! screening_lookup.clash(ms->beads())) {
            return ms;
        }
    }
    return ms;
}

void
DockMotifApp::_search() {
    auto ms = rotations_[0];
    auto rng = RandomNumberGenerator();
    int bound = 2;
    auto p = Point();

    float closest = 1000;
    float dist = 1000;
    for(auto const & b : ms->beads()) {
        dist = b.distance(center);
        if (dist < closest)  { closest = dist; }
    }

    auto current = closest;
    auto best = closest;
    auto next = closest;
    auto mc = MonteCarlo();
    auto best_ms = std::make_shared<MotifState>(*ms);
    auto current_pos = 0;

    for (int i = 0; i < 10000; i++) {
        p = get_random_point(rng, bound);
        ms->move(p);

        if(lookup_.clash(ms->beads())) {
            ms->move(-p);
        }

        closest = 1000;
        for(auto const & b : ms->beads()) {
            dist = b.distance(center);
            if (dist < closest)  { closest = dist; }
        }

        if(! mc.accept(current, dist)) {
            ms->move(-p);
        }

        current = dist;
        if(current < best) {
            best = current;
            best_ms = std::make_shared<MotifState>(*ms);
        }
    }

    auto m = RM::instance().get_motif_from_state(best_ms);
    std::cout << best << std::endl;
    m->to_pdb("docked.pdb");


}

Point
get_random_point(
        RandomNumberGenerator & rng,
        int bound) {
    auto x = bound - rng.rand()*2*bound;
    auto y = bound - rng.rand()*2*bound;
    auto z = bound - rng.rand()*2*bound;
    return Point(x, y, z);
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

Matrix
DockMotifApp::_rotation_about_x_axis(
        float degrees) {
    float a = 3.14/180*degrees;
    return Matrix(1.0, 0.0, 0.0,
                  0.0, cos(a), sin(a),
                  0.0, -sin(a), cos(a));
}

Matrix
DockMotifApp::_rotation_about_y_axis(
        float degrees) {
    float b = 3.14/180*degrees;
    return Matrix(cos(b), 0.0, -sin(b),
                  0.0, 1.0, 0.0,
                  sin(b), 0.0, cos(b));
}

Matrix
DockMotifApp::_rotation_about_z_axis(
        float degrees) {
    float g = 3.14/180*degrees;
    return Matrix(cos(g), sin(g), 0.0,
                  -sin(g), cos(g), 0.0,
                  0.0, 0.0, 1.0);
}

void
DockMotifApp::_precompute_rotations(
        MotifStateOP ms) {
    auto r_x = _rotation_about_x_axis(30);
    auto r_y = _rotation_about_y_axis(30);
    auto r_z = _rotation_about_z_axis(30);

    auto p = Point(0, 0, 0);

    ms->move(-ms->end_states()[0]->d());

    rotations_ = MotifStateOPs();

    auto count = 0;
    for(int i = 0; i < 12; i++) {
        ms->transform(r_x, p);
        for(int j = 0; j < 12; j++) {
            ms->transform(r_y, p);
            for(int k = 0; k < 12; k++) {
                ms->transform(r_z, p);
                rotations_.push_back(std::make_shared<MotifState>(*ms));
                //auto m = RM::instance().get_motif_from_state(ms);
                //m->to_pdb("frame."+std::to_string(count)+".pdb");
                //count++;

                //if(count > 10) { exit(0); }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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