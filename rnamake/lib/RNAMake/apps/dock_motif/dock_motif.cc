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

DockMotifApp::DockMotifApp() : base::Application() {}

void
DockMotifApp::setup_options() {

    add_option("scaffold", String(""), base::OptionType::STRING, true);
    add_option("motif", String(""), base::OptionType::STRING, true);
    add_option("origin", String(""), base::OptionType::STRING, false);
    add_option("ligand", String(""), base::OptionType::STRING, false);
}

void
DockMotifApp::parse_command_line(
        int argc,
        const char ** argv) {

    base::Application::parse_command_line(argc, argv);

}

void
DockMotifApp::run() {

    auto mf = motif::MotifFactory();
    auto added_helix = mf.added_helix();
    auto motif_path = get_string_option("motif");

    auto scaffold = resources::Manager::instance().get_structure(get_string_option("scaffold"), "scaffold");

    resources::Manager::instance().add_motif(get_string_option("motif"), "motif", util::MotifType::HAIRPIN);
    auto motif = resources::Manager::instance().motif("motif");
    motif->ends()[0]->flip();

    helix_ = resources::Manager::instance().motif_state("HELIX.IDEAL.2");

    auto start_ms = motif->get_state();
    _precompute_rotations(start_ms);

    auto screening_lookup = util::StericLookup(1.0, 5.0, 6);
    lookup_ = util::StericLookup(1.0, 4.0, 6);
    auto bead_centers = math::Points();
    for(auto const & b: scaffold->get_beads()) {
        bead_centers.push_back(b.center());
    }
    screening_lookup.add_points(bead_centers);
    lookup_.add_points(bead_centers);

    auto center = math::Point();
    if(get_string_option("ligand") != "") {
        auto ligand = _parse_ligand_for_center_coords();
        center_ = ligand->center();
    }

    results_ = MotifStateandScores();

    for (int i = 0; i < 1000; i++) {
        std::cout << i << std::endl;
        _get_starting_state(screening_lookup, center_);
        _search();

    }

    std::sort(results_.begin(), results_.end(), sort_by_score);

    int i = 0;
    for(auto const & result : results_) {
        std::cout << i << " : " << result.score << std::endl;
        auto m = resources::Manager::instance().get_motif_from_state(result.ms);
        m->to_pdb("docked." + std::to_string(i) + ".pdb");
        i++;
        if(i > 9) { break; }
    }

    //resources::Manager::instance().add_motif()


    //points_to_pdb("grid.pdb", valid_start_points);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

motif::MotifStateOP
DockMotifApp::_get_starting_state(
        util::StericLookup & screening_lookup,
        math::Point const & center) {

    auto ms = rotations_[0];
    auto rng = util::RandomNumberGenerator();
    auto p = math::Point();
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
    auto pos = 0;
    auto next_pos = 0;
    auto ms = rotations_[pos];
    auto last_ms = ms;
    auto rng = util::RandomNumberGenerator();
    int bound = 2;
    auto p = math::Point();

    auto current = _score(ms);
    auto best = current;
    auto next = current;
    auto mc = util::MonteCarlo();
    auto best_ms = std::make_shared<motif::motif::MotifState>(*ms);

    for (int i = 0; i < 100000; i++) {
        last_ms = ms;
        ms = rotations_[rng.randrange(rotations_.size()-1)];
        p = get_random_point(rng, bound);
        ms->move(p);

        if(lookup_.clash(ms->beads())) {
            ms = last_ms;
            continue;
        }

        next = _score(ms);
        if(! mc.accept(current, next)) {
            ms = last_ms;
            continue;
        }

        align_motif_state(ms->end_states()[0], helix_);
        if(lookup_.clash(helix_->beads())) {
            ms = last_ms;
            continue;
        }


        current = next;
        if(current < best) {
            best = current;
            best_ms = std::make_shared<motif::motif::MotifState>(*ms);
        }
    }

    results_.push_back(MotifStateandScore{best_ms, best});

    //auto m = resources::Manager::instance().get_motif_from_state(best_ms);
    //std::cout << best << std::endl;
    //m->to_pdb("docked.pdb");
}

float
DockMotifApp::_score(
        motif::MotifStateOP ms) {
    /*float avg = 0;
    for(auto const & b : ms->beads()) {
        avg += b.distance(center_);
    }
    return avg /= (float)ms->beads().size();*/

    float best = 1000;
    float dist = 0;
    for(auto const & b : ms->beads()) {
        dist = b.distance(center_);
        if(dist < best) { best = dist; }
    }
    return best;

}


math::Point
get_random_point(
        util::RandomNumberGenerator & rng,
        int bound) {
    auto x = bound - rng.rand()*2*bound;
    auto y = bound - rng.rand()*2*bound;
    auto z = bound - rng.rand()*2*bound;
    return math::Point(x, y, z);
}

math::Point
DockMotifApp::_calc_motif_center(
        motif::MotifOP m ) {
    auto m_center = math::Point();
    for(auto b : m->beads()) {
        m_center += b.center();
    }
    m_center /= m->beads().size();
    return m_center;
}

structure::ResidueOP
DockMotifApp::_parse_ligand_for_center_coords() {
    auto ligand_pdb_path = get_string_option("ligand");
    auto pdb_parser = structure::PDBParser();
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

math::Matrix
DockMotifApp::_rotation_about_x_axis(
        float degrees) {
    float a = 3.14/180*degrees;
    return math::Matrix(1.0, 0.0, 0.0,
                  0.0, cos(a), sin(a),
                  0.0, -sin(a), cos(a));
}

math::Matrix
DockMotifApp::_rotation_about_y_axis(
        float degrees) {
    float b = 3.14/180*degrees;
    return math::Matrix(cos(b), 0.0, -sin(b),
                  0.0, 1.0, 0.0,
                  sin(b), 0.0, cos(b));
}

math::Matrix
DockMotifApp::_rotation_about_z_axis(
        float degrees) {
    float g = 3.14/180*degrees;
    return math::Matrix(cos(g), sin(g), 0.0,
                  -sin(g), cos(g), 0.0,
                  0.0, 0.0, 1.0);
}

void
DockMotifApp::_precompute_rotations(
        motif::MotifStateOP ms) {
    auto r_x = _rotation_about_x_axis(30);
    auto r_y = _rotation_about_y_axis(30);
    auto r_z = _rotation_about_z_axis(30);

    auto p = math::Point(0, 0, 0);

    ms->move(-ms->end_states()[0]->d());

    rotations_ = motif::MotifStateOPs();

    auto count = 0;
    for(int i = 0; i < 12; i++) {
        ms->transform(r_x, p);
        for(int j = 0; j < 12; j++) {
            ms->transform(r_y, p);
            for(int k = 0; k < 12; k++) {
                ms->transform(r_z, p);
                rotations_.push_back(std::make_shared<motif::motif::MotifState>(*ms));
                //auto m = resources::Manager::instance().get_motif_from_state(ms);
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
    //String base_path = base::base_dir() + "/rnamake/lib/RNAMake/apps/simulate_tectos/resources/";

    auto app = DockMotifApp();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;

}