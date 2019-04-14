//
// Created by Joseph Yesselman on 3/22/19.
//

#ifndef RNAMAKE_NEW_SOLUTION_FILTER_H
#define RNAMAKE_NEW_SOLUTION_FILTER_H

#include <memory>

#include <base/types.h>
#include <base/string.h>

namespace motif_search {

class SolutionFilter {
public:
    SolutionFilter() {}

    virtual
    ~SolutionFilter() {}

    virtual
    SolutionFilter *
    clone() const = 0;

public:

    virtual
    bool
    accept(Strings const &) = 0;
};

typedef std::shared_ptr<SolutionFilter> SolutionFilterOP;

class NoExclusionFilter : public  SolutionFilter {
public:
    NoExclusionFilter(): SolutionFilter() {}

    ~NoExclusionFilter() {}

    SolutionFilter *
    clone() const {
        return new NoExclusionFilter(*this);
    }

public:

    bool
    accept(
            Strings const & motif_names) { return true; }
};

// removes helices in the same position of the same length mostly useful for flex helices
class RemoveDuplicateHelices : public  SolutionFilter {
public:
    RemoveDuplicateHelices(): SolutionFilter() {
        seen_ = StringIntMap();
    }

    ~RemoveDuplicateHelices() {}

    SolutionFilter *
    clone() const {
        return new RemoveDuplicateHelices(*this);
    }

public:

    bool
    accept(
            Strings const & motif_names) {
        key_ = "";
        for(auto const & m_name : motif_names) {
            if(m_name.substr(0, 10) == "HELIX.FLEX") {
                auto spl = base::split_str_by_delimiter(m_name, ".");
                key_ += spl[0] + "." + spl[1] + "." + spl[2] + ";";
            }
            else {
                key_ += m_name + ";";
            }
        }
        if(seen_.find(key_) != seen_.end()) { return false; }
        seen_[key_] = 1;

        return true;
    }

private:
    StringIntMap seen_;
    String key_;
};

}

#endif //RNAMAKE_NEW_SOLUTION_FILTER_H





























