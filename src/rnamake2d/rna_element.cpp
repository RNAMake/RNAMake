#include <rnamake2d/rna_element.h>
#include <iostream>

namespace rnamake2d {

    Ints
    get_pairmap_from_secstruct(String const & secstruct) {
       auto pair_stack = Ints{};
       auto pairs_array = Ints(secstruct.size(),-1);

       for(auto ii = 0; ii < secstruct.size(); ++ii) {
            if(secstruct[ii] == '(') {
                pair_stack.push_back(ii);
            } else if (secstruct[ii] == ')') {
                const auto index = *pair_stack.rbegin();
                pair_stack.pop_back();
                pairs_array[index] = ii;
                pairs_array[ii] = index;
            }
       }
       return pairs_array;
    }

    void
    get_rna_elements_from_secstruct_recursive(Ints const& pairs_array,
                                              int start_index,
                                              int end_index,
                                              RNAElems& elements,
                                              int last_pair_start,
                                              int last_pair_end,
                                              RNAElement* last_parent) {
        auto new_element = RNAElement{};
        new_element.type_ = RNAELEMENT::LOOP;

        auto ii = start_index;

        while(ii <= end_index) {
            if (pairs_array[ii] < 0) {
                new_element.indices_.push_back(ii);

                if (last_pair_start >= 0) {
                    auto stack_element = RNAElement();
                    stack_element.type_ = RNAELEMENT::STACK;
                    stack_element.parent_ = last_parent;
                    if (last_parent != nullptr) {
                        last_parent->children_.push_back(&stack_element);
                    }

                    for(auto jj = last_pair_start; jj < ii; ++jj) {
                        stack_element.indices_.push_back(jj);
                        stack_element.indices_.push_back(pairs_array[jj]);
                    }

                    elements.push_back(stack_element);
                    last_pair_start = -999;
                    last_pair_end = -999;
                    new_element.parent_ = &stack_element;
                    stack_element.children_.push_back(&new_element);
                }
                ++ii;
            }else if( ii < pairs_array[ii] ) {
                if( ii == start_index and pairs_array[ii] == end_index) {
                    if (last_pair_start < 0 ) {
                        last_pair_start = ii;
                        last_pair_end = pairs_array[ii];
                        last_parent = &new_element;
                    }
                    get_rna_elements_from_secstruct_recursive(
                            pairs_array,
                            ii + 1,
                            pairs_array[ii] - 1,
                            elements,
                            last_pair_start,
                            last_pair_end,
                            last_parent
                    );
                    return;
                } else {
                    if( last_pair_start >= 0 ) {
                        auto stack_element = RNAElement();
                        stack_element.type_ = RNAELEMENT::STACK;
                        stack_element.parent_ = last_parent;
                        if (last_parent != nullptr ) {
                            last_parent->children_.push_back(&stack_element);
                        }

                        for(auto jj = last_pair_start; jj < ii; ++jj) {
                            stack_element.indices_.push_back(jj);
                            stack_element.indices_.push_back(pairs_array[jj]);
                        }

                        elements.push_back(stack_element);
                        new_element.parent_ = &stack_element;
                        stack_element.children_.push_back(&new_element);
                        last_pair_start = -999;
                        last_pair_end = -999;
                    }

                    get_rna_elements_from_secstruct_recursive(
                            pairs_array,
                            ii + 1,
                            pairs_array[ii] - 1,
                            elements,
                            ii,
                            pairs_array[ii],
                            &new_element
                    );
                }
                ii = pairs_array[ii] + 1;
            } else {
                LOGE<<start_index<<" "<<end_index;
                LOGE<<"error : should not get here "<<ii<<" "<<pairs_array[ii];
                exit(1);
            }
        }

        elements.push_back(new_element);
    }

    RNAElems
    get_rna_elemnts_from_secstruct(String const& secstruct) {
       auto elements = RNAElems{};
       auto pairs_array = get_pairmap_from_secstruct(secstruct);

        get_rna_elements_from_secstruct_recursive(
                pairs_array, 0, secstruct.size() - 1, elements, -999,-999, nullptr
                );

       return elements;
    }

    Strings
    findall(String const& source, std::regex const& pattern)  {
        auto sm = std::smatch();
        std::regex_match(source, sm, pattern);
        auto matches = Strings{};
        for(auto ii = 0; ii < sm.size(); ++ii) {
            matches.push_back(sm[ii].str());
        }
        return matches;
    }




} // namespace rnamake2d