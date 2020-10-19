#include <rnamake2d/SSMotif.h>


namespace rnamake2d {

    void
    find_buffers(Motif2DOPs& motifs, int size) {
        std::vector<int> features(size);

        for(auto& motif : motifs) {
            if(motif->mtype() != util::MotifType::HELIX) {
                continue;
            }
            for(const auto& strand : motif->strands)  {
                for(const auto& pos : strand) {
                    features[pos] = motif->buffer_; //LOL
                }
            }
        }

        for(auto& motif : motifs) {
            auto best_buffer(std::numeric_limits<int>::max());
            for(const auto& strand : motif->strands) {
                best_buffer = std::min(best_buffer, features[*strand.cbegin()]);
                best_buffer = std::min(best_buffer, features[*strand.crbegin()]);
            }
            motif->buffer_ = best_buffer;
        }

    }


    Motif2DOPs
    parse_to_motif2ds(secondary_structure::PoseOP const & pose , String const& dot_bracket ) {
        auto motif2ds = Motif2DOPs{};
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

        for(const auto& motif : pose->motifs()) {
            auto strands = std::vector<Ints>();
            for (auto& res : motif->chains()) {
                strands.emplace_back();
                auto last = strands.rbegin();
                for(auto& bases : *res) {
                    last->push_back(bases->num()-1) ;
                }
            }
            switch( motif->mtype() ) {
                case util::MotifType::HAIRPIN: {
                    motif2ds.push_back(std::make_shared<rnamake2d::Hairpin>(strands));
                    continue;
                }
                case util::MotifType::SSTRAND: {
                    motif2ds.push_back(std::make_shared<rnamake2d::SingleStrand>(strands));
                    continue;
                }
                case util::MotifType::HELIX : {
                    continue;
                }
                case util::MotifType::NWAY : {
                    motif2ds.push_back(std::make_shared<rnamake2d::Junction>(strands));
                    continue;
                }
                case util::MotifType::TWOWAY: {
                    motif2ds.push_back(std::make_shared<rnamake2d::Junction>(strands));
                    continue;
                }
                default: {
                    throw std::runtime_error("Invalid motif code");
                }
            }
        }

        // edge case for single strands
        auto chain = Ints{};
        for(auto ii = 0; ii < dot_bracket.size(); ++ii) {
            const auto& ch = dot_bracket[ii];
            if (ch != '.') {
                if(!chain.empty() || ii == dot_bracket.size() - 1) {
                    bool used(false);
                    for(const auto& m : motif2ds) {
                        if(m->contains(ii-1)) {
                            used = true; break;
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

        find_buffers( motif2ds, dot_bracket.size() );

        return motif2ds;
    }
}

