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


private:
    MotifStateTreeOP mst_;
    CartesianProduct<Strings> pair_iterator_;
    std::vector<Ints> disallowed_num_sequences_;
    String seq1_, seq2_, seq_full_, structure_;
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
    MotifStateOPs
    get_motifs_from_seq_and_ss(
            String const &,
            String const &);

    void
    generate_helices(
            int);

    MotifOP
    get_avg_helix_new(
            int);

    MotifOP
    get_avg_helix(
            int);

    int
    convert_char_to_res_code(
            char);

    bool
    find_seq_violations(
            Ints const &);

    bool
    find_gc_strech(
            Ints const &);

private:
    std::vector<Ints> disallowed_num_sequences_;
};



#endif //TEST_BUILD_FLEX_HELICES_H































