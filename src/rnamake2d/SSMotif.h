#ifndef __SSMOTIF_H__
#define __SSMOTIF_H__

#include <unordered_set>
#include <memory>
#include <vector>
#include <string_view>
#include <regex>


#include <base/types.h>
#include <plog/Log.h>
#include <secondary_structure/secondary_structure_parser.h>

namespace rnamake2d {
    struct Motif {
        // what are motifs going to be important for? knowing what sequence they have and the helix size
        std::string_view full_sequence_; // this is owned by the design
        std::string_view full_mutated_sequence_;
        std::map<int,int> pairs;
        String sequence; // this is owned by the Motif
        String dot_bracket; // this is owned by the Motif
        String token_;
        int last_change;
        int size_;
        std::vector<Ints> strands;
        std::unordered_set<int> nts_;
        int buffer_;
        util::MotifType mtype_;

        private:
        void
        build_sequence_(bool mutant=false)  {
            if(!sequence.empty()) {
                sequence.clear();
            }
            for(const auto& span : strands) {
                for(const auto index : span)  {
                    if(mutant)  {
                        sequence += full_mutated_sequence_[index];
                    } else {
                        sequence += full_sequence_[index];
                    }
                }
                sequence += '&';
            }
            sequence.pop_back();
        }

    public:
        String
        token() const {
            return token_;
        }

        void
        update(bool mutant=false)  {

        }

        virtual
        void
        full_sequence(String const& full_sequence )  {
            full_sequence_ = full_sequence;
            build_sequence_();
        }

        virtual
        ~Motif() = default;

        void
        show() {
            std::cout<<token_<<' ';
            for(auto& strand :strands ) {
                for(auto& n : strand) {
                    std::cout<<n<<",";
                }
            }
            std::cout<<" Buffer: "<<buffer_<<", Seq: ";
            if(!sequence.empty()) {
                std::cout<<sequence;
            }
            std::cout<<"\n";
        }

        int
        buffer() const {
            return buffer_;
        }

        bool
        contains(int num) {
            return nts_.find(num) != nts_.end();
        }

        util::MotifType
        mtype() const {
            return mtype_;
        }

    protected:
        void
        setup_nts_() {
            for(const auto& strand : strands) {
                for(const auto& pos : strand) {
                    nts_.insert(pos);
                }
            }
        }
    };

    struct Hairpin : Motif {
        int size_ = -1;
        explicit Hairpin(std::vector<Ints>& strands) {
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
        }
    };

    struct Junction : Motif {
        int num_branches_ = 0;
        Ints gap_sizes;
        int gc{0}, au{0}, gu{0}, unknown{0};
        explicit Junction(std::vector<Ints>& strands) {
            if(strands.size() < 2) {
                std::cout<<"Junction must have at LEAST 2 strands\n";
                exit(1);
            }
            this->strands = std::move(strands);
            token_ = "Junction" + std::to_string(this->strands.size());
            mtype_ = util::MotifType::NWAY;
            for(const auto& strand : this->strands) {
                token_ += '_';
                token_ += std::to_string(*strand.crbegin() - *strand.cbegin());
            }
            setup_nts_();
        }

        void
        full_sequence(String const & sequence)  override {
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

    };

    struct Helix : Motif {
        explicit Helix(std::vector<Ints>& strands) {
            if(strands.size() != 2) {
                std::cout<<"Helix must be 2 strand\n";
                exit(1);

            }
            mtype_ = util::MotifType::HELIX;
            this->strands = std::move(strands);
            buffer_ = this->strands[0].size();
            size_ = this->strands[0].size();
            token_ = "Helix" + std::to_string(size_);
            setup_nts_();
        }
    };

    struct SingleStrand : Motif {

        explicit SingleStrand(std::vector<Ints>& strands) {
            if(strands.size() != 1) {
                std::cout<<"SingleStrand must be 1 strand\n";
                exit(1);
            }
            this->strands = std::move(strands);
            size_ = this->strands.size();
            token_ = "SingleStrand" + std::to_string(size_);
            setup_nts_();
        }
    };
}

using Motif2D = rnamake2d::Motif;
using Motif2Ds = std::vector<Motif2D>;

using Motif2DOP = std::shared_ptr<Motif2D>;
using Motif2DOPs = std::vector<Motif2DOP>;

namespace  rnamake2d {

    Motif2DOPs
    parse_to_motif2ds( secondary_structure::PoseOP const&, String const& );

    void
    find_buffers( Motif2DOPs&, int );
}

#endif // __SSMOTIF_H__