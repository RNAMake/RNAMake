//
// Created by Joseph Yesselman on 3/3/19.
//

#ifndef TEST_BUILD_FLEX_HELICES_H
#define TEST_BUILD_FLEX_HELICES_H

#include <stdio.h>
#include "base/application.hpp"
#include "util/cartesian_product.h"
#include "secondary_structure/secondary_structure_parser.h"
#include "motif/motif.h"
#include "resources/resource_manager.h"
#include <motif_data_structures/motif_state_tree.h>


class HelixStructureIterator {
public:
    HelixStructureIterator() {
        // disallow all 4 residue of the same type in an row "AAAA" etc
        disallowed_num_sequences_ = std::vector<Ints>();
        for(int i = 0; i < 4; i++) {
            auto disallowed_sequence = Ints(4);
            for(int j = 0; j < 4; j++) { disallowed_sequence[j] = i; }
            disallowed_num_sequences_.push_back(disallowed_sequence);
        }

        mst_ = std::make_shared<MotifStateTree>();
        mst_->set_option_value("sterics", false);
    }

    ~HelixStructureIterator() {}

public:
    void
    setup(
            int length) {
        auto pairs = std::vector<Strings>();
        pairs.push_back(Strings{"A", "U"});
        pairs.push_back(Strings{"U", "A"});
        pairs.push_back(Strings{"C", "G"});
        pairs.push_back(Strings{"G", "C"});

        auto all_pairs = std::vector<std::vector<Strings>>((int)(length));
        for(int i = 0; i < all_pairs.size(); i++) {
            all_pairs[i] = pairs;
        }
        pair_iterator_ = CartesianProduct<Strings>(all_pairs);
        motifs_ = MotifStateOPs(length-1);
        num_seq_ = Ints(length);

        for(int i = 0; i < length; i++) { structure_ += "("; }
        structure_ += "&";
        for(int i = 0; i < length; i++) { structure_ += ")"; }

    }

    MotifStateTreeOP
    next() {
        auto done = false;
        while(!done && !pair_iterator_.end()) {
            auto & current = pair_iterator_.next();
            seq1_ = ""; seq2_ = "";
            for (auto const & p : current) {
                seq1_ += p[0];
                seq2_ = p[1] + seq2_;
            }

            for(int i = 0 ; i < num_seq_.size(); i++) {
                num_seq_[i] = convert_char_to_res_code(seq1_[i]);
            }
            if(find_seq_violations(num_seq_) || find_gc_strech(num_seq_)) { continue; }

            seq_full_ = seq1_ + "&" + seq2_;
            get_motifs_from_seq_and_ss(seq_full_, structure_, motifs_);

            mst_->remove_node_level(0);
            for(auto const & m : motifs_) {
                mst_->add_state(m);
            }

            break;
        }

        return mst_;
    }

    bool
    end() {
        return pair_iterator_.end();
    }

public:
    String const &
    get_current_sequence() {
        return seq_full_;
    }

    MotifStateTreeOP
    get_tree_with_sequence(
            String const & seq) {
        if(seq.size() != structure_.size()) {
            throw std::runtime_error("cannot build new tree seq and structure not the same size");
        }
        get_motifs_from_seq_and_ss(seq, structure_, motifs_);

        auto new_mst = std::make_shared<MotifStateTree>();
        for(auto const & m : motifs_) {
            new_mst->add_state(m);
        }
        return new_mst;
    }

private:

    void
    get_motifs_from_seq_and_ss(
            String const & seq,
            String const & ss,
            MotifStateOPs & motifs) {
        auto parser = sstruct::SecondaryStructureParser();
        auto ss_motifs = parser.parse_to_motifs(seq, ss);

        auto start = 0;
        auto motif = MotifStateOP(nullptr);
        int i = 0;
        for(auto const & m : ss_motifs) {
            //basepair step
            if(m->mtype() == MotifType::HELIX) {
                motif = RM::instance().bp_step(m->end_ids()[0])->get_state();
                motifs[i] = motif;
            }
            else {
                throw std::runtime_error("only helices are allowed");
            }
            i++;
        }
    }

    int
    convert_char_to_res_code(
            char c) {
        if     (c == 'A') { return 0; }
        else if(c == 'C') { return 1; }
        else if(c == 'G') { return 2; }
        else if(c == 'U') { return 3; }
        else if(c == 'T') { return 3; }
        else if(c == 'N') { return -1; }
        else {
            throw std::runtime_error("incorrect character for secondary string");
        }
        return -1;
    }

    bool
    find_seq_violations(
            Ints const & num_seq) {
        auto pos = 0;
        for (int i = 0; i < num_seq.size(); i++) {
            for (int j = 0; j < disallowed_num_sequences_.size(); j++) {
                if (i < disallowed_num_sequences_[j].size()-1) { continue; }
                auto match = true;
                pos = i - (disallowed_num_sequences_[j].size()-1);
                if(pos < 0) { continue; }
                for (auto const & e : disallowed_num_sequences_[j]) {
                    if (num_seq[pos] != e) {
                        match = false;
                        break;
                    }
                    pos += 1;
                }
                if (match) {
                    return true;
                }
            }
        }
        return false;
    }

    bool
    find_gc_strech(
            Ints const & num_seq) {
        int count = 0;
        for(auto const & r : num_seq) {
            if(r == 1 || r == 2)  { count += 1; }
            else                  { count = 0;  }
            if(count > 2) { return true; }
        }
        return false;
    }


private:
    MotifStateTreeOP mst_;
    CartesianProduct<Strings> pair_iterator_;
    std::vector<Ints> disallowed_num_sequences_;
    String seq1_, seq2_, seq_full_, structure_;
    Ints num_seq_;
    MotifStateOPs motifs_;
};


class BuildFlexHelicesAppException : public std::runtime_error {
public:
    BuildFlexHelicesAppException(
            String const & message):
            std::runtime_error(message)
    {}
};

class BuildFlexHelicesApp : public Application {
public:
    BuildFlexHelicesApp();

    ~BuildFlexHelicesApp() {}

public:
    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    void
    generate_helices(
            int);

    MotifOP
    get_avg_helix(
            int);

private:
    MotifFactory mf_;
};



#endif //TEST_BUILD_FLEX_HELICES_H































