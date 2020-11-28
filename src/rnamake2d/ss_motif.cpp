#include <rnamake2d/ss_motif.h>


namespace rnamake2d {
    void
    Motif2D::build_sequence_()  {
        if(!sequence.empty()) {
            sequence.clear();
        }

        for(const auto& span : strands) {
            for(const auto index : span)  {
                sequence += full_sequence_[index];
            }
            sequence += '&';
        }

        sequence.pop_back();
    }

    bool
    Motif2D::operator==(const Motif& other) const {

        return mtype_ == other.mtype_               &&
               size_ == other.size_                 &&
               buffer_ == other.buffer_             &&
               sequence == other.sequence           &&
               dot_bracket == other.dot_bracket     &&
               comp_vectors(strands, other.strands) &&
               comp_usets(nts_, other.nts_)         &&
               token_ == other.token_;

    }

    void
    Motif2D::setup_nts_() {
        for(const auto& strand : strands) {
            for(const auto& pos : strand) {
                nts_.insert(pos);
            }
        }
    }

    Hairpin::Hairpin(std::vector<Ints>& strands, String const& sequence) {
        if(strands.size() != 1) {
            std::cout<<"Hairpin must be 1 strand\n";
            exit(1);
        }
        this->strands = std::move(strands);

        size_ = this->strands[0].size() - 2;
        token_ = "Hairpin" + std::to_string(size_);
        buffer_ = 1;
        mtype_ = util::MotifType::HAIRPIN;
        setup_nts_();

        if(!sequence.empty()) {
            full_sequence(sequence);
        }
    }

    Junction::Junction(std::vector<Ints> &strands, String const& sequence) {
        if(strands.size() < 2) {
            std::cout<<"Junction must have at LEAST 2 strands\n";
            exit(1);
        }
        this->strands = std::move(strands);
        token_ = "Junction" + std::to_string(this->strands.size());
        mtype_ = util::MotifType::NWAY;
        for(const auto& strand : this->strands) {
            token_ += '_';
            token_ += std::to_string(*strand.crbegin() - *strand.cbegin() - 1);
        }
        setup_nts_();
        //std::cout<<sequence<<std::endl; exit(1);
        if(!sequence.empty()) {
            full_sequence(sequence);
        }
    }

    void
    Junction::full_sequence(const String &sequence) {
        this->Motif::full_sequence(sequence);
        auto pairs = std::vector<std::pair<int,int>>{{*strands.begin()->begin(),*strands.rbegin()->rbegin()}};

        auto it = sequence.find('&');
        while(it != String::npos)  {
            pairs.emplace_back(it-1,it+1);
            it = sequence.find('&',it+1);
        }

        for(const auto& pr : pairs) {
            String base_pair;
            base_pair += sequence[pr.first] + sequence[pr.second];
            if(base_pair == "AU" || base_pair == "UA") {
                ++au;
            } else if (base_pair == "GC" || base_pair == "CG") {
                ++gc;
            } else if (base_pair == "GU" || base_pair == "UG") {
                ++gu;
            } else {
                ++unknown;
            }
        }
    }

    Helix::Helix(std::vector<Ints> &strands, String const& sequence) {
        if(strands.size() != 2) {
            LOGE<<"Helix must be 2 strand\n";
            exit(1);

        }
        mtype_ = util::MotifType::HELIX;
        this->strands = std::move(strands);
        buffer_ = this->strands[0].size();
        size_ = buffer_;
        token_ = "Helix" + std::to_string(size_);
        setup_nts_();

        if(!sequence.empty()) {
            full_sequence(sequence);
        }
    }

    SingleStrand::SingleStrand(std::vector<Ints> &strands, String const& sequence) {
        mtype_ = util::MotifType::SSTRAND;
        buffer_ = 0;
        if(strands.size() != 1) {
            std::cout<<"SingleStrand must be 1 strand\n";
            exit(1);
        }
        this->strands = std::move(strands);
        size_ = this->strands.size();
        token_ = "SingleStrand" + std::to_string(size_);
        setup_nts_();
        if(!sequence.empty()) {
            full_sequence(sequence);
        }
    }

    void
    find_buffers(HelixOPs& helices, JunctionOPs& junctions, SingleStrandOPs& singlestrands, HairpinOPs& hairpins, int size) {
        std::vector<int> features(size);
        for(auto& motif : helices) {
            for(const auto& strand : motif->strands)  {
                for(const auto& pos : strand) {
                    features[pos] = motif->buffer_; //LOL
                }
            }
        }

        for(auto& motif : junctions) {
            auto best_buffer(std::numeric_limits<int>::max());
            for(const auto& strand : motif->strands) {
                best_buffer = std::min(best_buffer, features[*strand.cbegin()]);
                best_buffer = std::min(best_buffer, features[*strand.crbegin()]);
            }
            motif->buffer_ = best_buffer;
        }

        for(auto& motif : singlestrands) {
            auto best_buffer(std::numeric_limits<int>::max());
            for(const auto& strand : motif->strands) {
                best_buffer = std::min(best_buffer, features.at(*strand.cbegin()));
                best_buffer = std::min(best_buffer, features.at(*strand.crbegin()));
            }
            motif->buffer_ = best_buffer;
        }

        for(auto& motif : hairpins) {
            auto best_buffer(std::numeric_limits<int>::max());
            for(const auto& strand : motif->strands) {
                best_buffer = std::min(best_buffer, features[*strand.cbegin()]);
                best_buffer = std::min(best_buffer, features[*strand.crbegin()]);
            }
            motif->buffer_ = best_buffer;
        }


    }

    static
    void setup_helices(secondary_structure::PoseOP const & pose, HelixOPs & motif2ds ) {
        for(const auto& motif : pose->helices()) {
            // helix will always have two strands
            if( motif->chains().size() != 2 ) {
                throw std::runtime_error("Helix size MUST be 2");
            }
            auto strands = std::vector<Ints>();
            for (auto& res : motif->chains()) {
                strands.emplace_back();
                auto last = strands.rbegin();
                for(auto& bases : *res) {
                    last->push_back(bases->num()-1) ;
                }
            }
            motif2ds.push_back(std::make_shared<rnamake2d::Helix>(strands));
        }
    }

    static
    void
    setup_hairpins_singlestrands_junctions(secondary_structure::PoseOP const & pose, HairpinOPs& hairpins, SingleStrandOPs& singlestrands, JunctionOPs& junctions ) {

        for(const auto& motif : pose->motifs()) {
            auto strands = std::vector<Ints>();
            for (auto& res : motif->chains()) {
                auto chain = Ints{};
                for(auto& bases : *res) {
                    chain.push_back(bases->num() - 1);
                }
                strands.push_back(std::move(chain));
            }
            const auto m_type = motif->mtype();
            if(m_type == util::MotifType::HAIRPIN) {
                hairpins.push_back(std::make_shared<rnamake2d::Hairpin>(strands));
            } else if (m_type == util::MotifType::SSTRAND) {
                singlestrands.push_back(std::make_shared<rnamake2d::SingleStrand>(strands));
            } else if (m_type == util::MotifType::NWAY || m_type == util::MotifType::TWOWAY) {
                junctions.push_back(std::make_shared<rnamake2d::Junction>(strands));
            }

        }
    }

    static
    void
    setup_edge_case_singlestrands(secondary_structure::PoseOP const & pose, HairpinOPs const& hairpins, HelixOPs const& helices, JunctionOPs const& junctions, SingleStrandOPs & motif2ds, const String& dot_bracket) {
        auto chain = Ints{};
        for(auto ii = 0; ii < dot_bracket.size(); ++ii) {
            const auto& ch = dot_bracket[ii];
            if (ch != '.') {
                if(!chain.empty() || ii == dot_bracket.size() - 1) {
                    bool used(false);
                    for(const auto& m : hairpins) {
                        if(m->contains(ii-1)) {
                            used = true;
                            break;
                        }
                    }
                    for(const auto& m : helices) {
                        if(m->contains(ii-1)) {
                            used = true;
                            break;
                        }
                    }
                    for(const auto& m : junctions) {
                        if(m->contains(ii-1)) {
                            used = true;
                            break;
                        }
                    }

                    if(!used) {
                        auto formatted_chain = std::vector<Ints>{chain};
                        motif2ds.push_back(std::make_shared<rnamake2d::SingleStrand>(formatted_chain));
                    }
                }
                chain.clear();
            } else {
                chain.push_back(ii);
            }
        }
    }

    std::tuple<HairpinOPs ,HelixOPs,  JunctionOPs, SingleStrandOPs>
    parse_to_motif2ds(secondary_structure::PoseOP const & pose , String const& dot_bracket, String const& sequence ) {
        //auto motif2ds = Motif2DOPs{};

        auto helices = HelixOPs{};
        auto hairpins = HairpinOPs{};
        auto junctions = JunctionOPs{};
        auto singlestrands = SingleStrandOPs{};
        setup_helices(pose, helices);

        setup_hairpins_singlestrands_junctions(pose, hairpins, singlestrands, junctions);

        setup_edge_case_singlestrands(pose, hairpins, helices, junctions, singlestrands, dot_bracket);

        find_buffers(helices, junctions, singlestrands, hairpins, dot_bracket.size());

        for(auto& hel : helices) {hel->full_sequence(sequence);}
        for(auto& hps : hairpins) {hps->full_sequence(sequence);}
        for(auto& junc : junctions) {junc->full_sequence(sequence);}
        for(auto& ss : singlestrands) {ss->full_sequence(sequence);}

        return std::make_tuple(std::move(hairpins), std::move(helices), std::move(junctions), std::move(singlestrands));
    }
}

