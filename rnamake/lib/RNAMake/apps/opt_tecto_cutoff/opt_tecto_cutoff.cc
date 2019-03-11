//
// Created by Joseph Yesselman on 2/19/18.
//
#include <algorithm>

#include "base/backtrace.hpp"
#include "math/hashing.h"
#include "math/stats.h"
#include "util/monte_carlo.h"
#include "util/random_number_generator.h"
#include "base/file_io.h"
#include "opt_tecto_cutoff/opt_tecto_cutoff.h"


OptTectoCutoff::OptTectoCutoff() : Application() {}

// application setups functions ////////////////////////////////////////////////////////////////////

void
OptTectoCutoff::setup_options() {
    //add_option("pdb", "", OptionType::STRING, true);
    add_option("histos", "", OptionType::STRING, true);
    add_option("data", "", OptionType::STRING, true);
    add_option("divide_set", false, OptionType::BOOL, false);
    add_option("divide_set_num", 2, OptionType::INT, false);
    add_option("constraints", "", OptionType::STRING, false);
    add_option("constraint_file", "", OptionType::STRING, false);


}

void
OptTectoCutoff::parse_command_line(
        int argc,
        const char ** argv) {
    Application::parse_command_line(argc, argv);
}


void
OptTectoCutoff::_setup() {

    auto lines =base::get_lines_from_file(get_string_option("data"));
    exp_dgs_ = std::vector<double>();
    constructs_ = Constructs();
    for(int i = 1; i < lines.size(); i++) {
        auto spl = base::split_str_by_delimiter(lines[i], ",");
        if(spl.size() < 3) { break; }
        constructs_.push_back(Construct{spl[0], spl[1], spl[2], spl[3], std::stod(spl[4]), i-1});
        exp_dgs_.push_back(std::stod(spl[4]));
        //if(exp_dgs_.size() > 100) { break; }
    }

    std::ifstream in;
    in.open(get_string_option("histos"), std::ios::binary);
    histos_ =  std::vector<SixDHistogram>();
    while (in.good()) {
        histos_.push_back(SixDHistogram(in));
        if(histos_.back().size() == 0) {
            exp_dgs_.erase(exp_dgs_.begin()+histos_.size()-1);
            constructs_.erase(constructs_.begin()+histos_.size()-1);
            histos_.pop_back();
        }
        if(histos_.size() == constructs_.size()) { break; }
    }
    in.close();
    /*for(int i = 0; i < histos_.size(); i++) {
        std::cout << i << " " << histos_[i].total_count() << " " << exp_dgs_[i] << std::endl;
    }

    std::cout << histos_.size() << " " << exp_dgs_.size() << std::endl;
    exit(0);*/

    avg_hit_counts_ = std::vector<double>(histos_.size());
    pred_dgs_ = std::vector<double>(histos_.size());
    norm_exp_dgs_ = std::vector<double>(histos_.size());

}

void
OptTectoCutoff::_divide_dataset() {
    std::sort(constructs_.begin(), constructs_.end(), CompareConstructDG());

    int step = get_int_option("divide_set_num");

    for(int j = 0; j < step; j++) {
        std::ofstream c_out ("constructs_"+std::to_string(j) + ".csv");
        std::ofstream h_out ("histos_"+std::to_string(j) + ".out", std::ios::binary);
        c_out << "fseq,fss,cseq,css,dG\n";
        int i = 0;
        for(i = 0; i < constructs_.size(); i+=step) {
            if(histos_[constructs_[i+j].pos].size() == 0) { continue; }
            c_out << constructs_[i+j].fseq << "," << constructs_[i+j].fss << ",";
            c_out << constructs_[i+j].cseq << "," << constructs_[i+j].css << ",";
            c_out << constructs_[i+j].dg << std::endl;
            histos_[constructs_[i+j].pos].output_binary(h_out, -1);
        }
        if(i+j < constructs_.size()) {
            c_out << constructs_[i+j].fseq << "," << constructs_[i+j].fss << ",";
            c_out << constructs_[i+j].cseq << "," << constructs_[i+j].css << ",";
            c_out << constructs_[i+j].dg << std::endl;
            histos_[constructs_[i+j].pos].output_binary(h_out, -1);
        }
        c_out.close();
        h_out.close();
    }

}

void
OptTectoCutoff::_get_scored_dataset() {
    auto constraints_str = get_string_option("constraints");

    auto spl = base::split_str_by_delimiter(constraints_str, ";");
    for (auto const & s : spl) {
        auto spl2 = base::split_str_by_delimiter(s, ",");
        auto pos = _parse_constraint_position(spl2[0]);
        auto lower = std::stod(spl2[1]);
        auto upper = std::stod(spl2[2]);
        constraints_[pos] = Real2{lower, upper};
    }

    auto score = _score(constraints_);

    std::cout << score << " " << r_ << " " << avg_diff_ << std::endl;

    std::ofstream out;
    out.open("scored.csv");
    out << "fseq,fss,cseq,css,dG,dG_normalized,dG_predicted,avg_hit_count\n";
    int i = 0;
    for (auto const & c : constructs_) {
        out << c.fseq << "," << c.fss << "," << c.cseq << "," << c.css << ",";
        out << c.dg << "," << norm_exp_dgs_[i] << "," << pred_dgs_[i] << "," << avg_hit_counts_[i] << std::endl;
        i++;
    }
    out.close();

}


int
OptTectoCutoff::_parse_constraint_position(
        String const & constraint_name) {
    if     (constraint_name == "x") { return 0; }
    else if(constraint_name == "y") { return 1; }
    else if(constraint_name == "z") { return 2; }
    else if(constraint_name == "a") { return 3; }
    else if(constraint_name == "b") { return 4; }
    else if(constraint_name == "g") { return 5; }
    else {
        throw std::runtime_error("unknown constraint name: " + constraint_name);
    }

}



void
OptTectoCutoff::_score_constraint_file() {
    auto constraint_file = get_string_option("constraint_file");
    std::ifstream in (constraint_file);
    auto line = String();
    int i = 0;

    std::ofstream out_sum ("constraint_scores.csv");
    out_sum << "x_min,x_max,y_min,y_max,z_min,z_max,a_min,a_max,b_min,b_max,g_min,g_max,r,avg_diff,lowest\n";

    while (in.good()) {
        getline(in, line);
        if(line[0] == 'x') {continue; }
        auto spl = base::split_str_by_delimiter(line, ",");
        int j = 0;
        for(int k = 0; k < 12; k += 2) {
            constraints_[j] = Real2{std::stod(spl[k]), std::stod(spl[k+1])};
            j += 1;
        }

        auto score = _score(constraints_);
        if(score > 100) { continue; }
        std::cout << i << " " << score << " " << r_ << " " << avg_diff_ << std::endl;

        auto lowest = 100000;
        for(auto const & hits : avg_hit_counts_) {
            if(hits < lowest) { lowest = hits; }
        }

        for(int k = 0; k < spl.size(); k++) {
            out_sum << spl[k] << ",";
        }

        out_sum <<  r_ << "," << avg_diff_ << "," << lowest << std::endl;
        //if(i > 100) { break; }
        i++;

    }
    out_sum.close();


}




double
OptTectoCutoff::_score(
        std::array<Real2, 6> const & constraints) {

    lowest_ = 1000000;
    for(int i = 0; i < histos_.size(); i++) {
        avg_hit_counts_[i] = histos_[i].within_constraints(constraints);
        if(avg_hit_counts_[i] < 1) { return 10000; }
        if(avg_hit_counts_[i] < lowest_) {
            lowest_ = avg_hit_counts_[i];
        }
    }
    auto avg_hit_count_mean = mean(avg_hit_counts_);


    double closest = 100000;
    int closest_i = 0;
    double diff = 0;
    for(int i = 0; i < histos_.size(); i++) {
        diff = abs(avg_hit_counts_[i] - avg_hit_count_mean);
        if(diff < closest) {
            closest = diff;
            closest_i = i;
        }
    }

    for(int i = 0; i < histos_.size(); i++) {
        pred_dgs_[i] = 1.9872041e-3*298*log(avg_hit_counts_[closest_i] / avg_hit_counts_[i]);
        norm_exp_dgs_[i] = exp_dgs_[i] - exp_dgs_[closest_i];
    }

    r_ = pearson_coeff(pred_dgs_, norm_exp_dgs_);
    avg_diff_ =  avg_unsigned_diff(pred_dgs_, norm_exp_dgs_);
    return  avg_diff_*4 + 10 / lowest_;

}

void
OptTectoCutoff::_vary_constraints(
        std::array<Real2, 6> const & constraints,
        std::array<Real2, 6> & new_constraints,
        RandomNumberGenerator & rng) {

    new_constraints = constraints;
    double val = 0;
    bool done = false;

    while(! done) {
        auto pos = rng.randrange(6);
        auto i = rng.randrange(2);
        auto diff = 1 + rng.randrange(3);
        auto dir = rng.randrange(2);

        if (pos < 3) { val = 0.25 * diff; }
        else { val = 5.0 * diff; }

        if (dir == 1) { val = -val; }

        if (pos < 3) {
            if(new_constraints[pos][i] + val > 6 || new_constraints[pos][i] + val < -6) { continue; }
        }
        else {
            if(new_constraints[pos][i] + val > 360 || new_constraints[pos][i] + val < 0) { continue; }
        }

        new_constraints[pos][i] += val;
        break;
    }


}


void
OptTectoCutoff::_get_initial_constraints(
        RandomNumberGenerator & rng) {
    constraints_ = std::array<Real2, 6>();

    bool done = false;
    while(!done) {
        for (int i = 0; i < 6; i++) {
            if (i < 3) {
                auto min = -rng.randrange(25) * .25;
                auto max = rng.randrange(25) * .25;
                constraints_[i][0] = min;
                constraints_[i][1] = max;
            }
            else {
                auto min = -rng.randrange(36) * 5.0 + 180;
                auto max = rng.randrange(36) * 5.0 + 180;
                constraints_[i][0] = min;
                constraints_[i][1] = max;
            }
        }
        auto score = _score(constraints_);
        for(auto const & c : constraints_) { std::cout << c[0] << " " << c[1] << " "; }
        std::cout << std::endl;
        if(score != 10000) { return; }
    }

}

void
OptTectoCutoff::run() {
    _setup();
    if(get_bool_option("divide_set")) {
        _divide_dataset();
        return;
    }
    if(get_string_option("constraints") != "") {
        _get_scored_dataset();
        return;
    }
    if(get_string_option("constraint_file") != "") {
        _score_constraint_file();
        return;
    }

    auto rng = RandomNumberGenerator();
    auto mc = MonteCarlo(0.05);
    _get_initial_constraints(rng);

    //constraints_ = {Real2{-5, 5}, Real2{-5, 5}, Real2{-5, 5}, Real2{10, 340}, Real2{130, 190}, Real2{170, 190}};
    //constraints_ = {Real2{-4.5, 4.0}, Real2{-4.5, 4.5}, Real2{-4.5, 4.5}, Real2{10, 340}, Real2{150, 190}, Real2{180, 200}};
    //constraints_ = {Real2{-4.5, 4.0}, Real2{-4.5, 4.5}, Real2{-4.5, 4.5}, Real2{10, 340}, Real2{150, 190}, Real2{140, 200}};

    new_constraints_ = std::array<Real2, 6>();


    auto current_score = _score(constraints_);
    auto new_score = current_score;

    std::cout << current_score << " " << r_ << " " << avg_diff_ << " " << lowest_ << std::endl;

    double best_score = 10000;
    auto best_constraints = constraints_;
    double best_avg_diff;
    double best_r;

    int heat_up = 0;
    for(int i = 0; i < 10000; i++) {
        if( i % 100 == 0) {
            constraints_ = best_constraints;
            current_score = best_score;
            std::cout << "STEPS: " << i << std::endl;
            heat_up = 0;
            mc.set_temperature(0.5);
        }

        heat_up += 1;
        if(heat_up > 20) {
            mc.set_temperature(0.05);
        }

        _vary_constraints(constraints_, new_constraints_, rng);
        new_score = _score(new_constraints_);
        if(mc.accept(current_score, new_score)) {
            constraints_ = new_constraints_;
            current_score = new_score;
            if(best_score > current_score) {
                best_score = current_score;
                best_constraints = constraints_;
                best_r = r_;
                best_avg_diff = avg_diff_;
                std::cout << "best: " << best_score << " " << best_r << " " << best_avg_diff << " " << lowest_ << std::endl;
                for(auto const & c : best_constraints) {
                    std::cout << c[0] << " " << c[1] << " ";
                }
                std::cout << std::endl;
            }
        }



    }




}


int main(int argc, const char * argv[]) {
    //must add this for all apps!
    std::set_terminate(print_backtrace);

    auto app = OptTectoCutoff();
    app.setup_options();
    app.parse_command_line(argc, argv);
    app.run();

    return 0;


}