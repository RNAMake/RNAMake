#ifndef __SSMOTIF_H__
#define __SSMOTIF_H__

#include <unordered_set>
#include <memory>
#include <vector>
#include <string_view>
#include <regex>
#include <tuple>

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
        build_sequence_(bool mutant=false);

    public:
        virtual
        bool
        is_helix() const {
            return false;
        }

    public:
        virtual
        bool
        is_hairpin() const {
            return false;
        }

    public:
        virtual
        bool
        is_junction() const {
            return false;
        }

    public:
        virtual
        bool
        is_singlestrand() const {
            return false;
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
        setup_nts_();

    };

    struct Hairpin : Motif {
        int size_ = -1;
        explicit Hairpin(std::vector<Ints>& strands);

    public:
        bool
        is_hairpin() const final {
            return true;
        }
    };

    struct Junction : Motif {
        int num_branches_ = 0;
        Ints gap_sizes;
        int gc{0}, au{0}, gu{0}, unknown{0};
        explicit Junction(std::vector<Ints>& strands);

        void
        full_sequence(String const & sequence)  override;

    public:
        bool
        is_junction() const final {
            return true;
        }
    };

    struct Helix : Motif {
        explicit Helix(std::vector<Ints>& strands);

    public:
        bool
        is_helix() const final {
            return true;
        }

    };

    struct SingleStrand : Motif {
        explicit SingleStrand(std::vector<Ints>& strands);

        public:
            bool
            is_singlestrand() const final {
                return true;
            }

    };
}

using Motif2D = rnamake2d::Motif;
using Motif2Ds = std::vector<Motif2D>;

using Motif2DOP = std::shared_ptr<Motif2D>;
using Motif2DOPs = std::vector<Motif2DOP>;

using HairpinOPs = std::vector<std::shared_ptr<rnamake2d::Hairpin>>;
using JunctionOPs = std::vector<std::shared_ptr<rnamake2d::Junction>>;
using HelixOPs = std::vector<std::shared_ptr<rnamake2d::Helix>>;
using SingleStrandOPs = std::vector<std::shared_ptr<rnamake2d::SingleStrand>>;

namespace  rnamake2d {

    std::tuple< HairpinOPs , HelixOPs, JunctionOPs, SingleStrandOPs>
    parse_to_motif2ds( secondary_structure::PoseOP const&, String const& );

    void
    find_buffers(HelixOPs& helices, JunctionOPs& junctions, SingleStrandOPs& singlestrands, HairpinOPs& hairpins, int size);
}

#endif // __SSMOTIF_H__