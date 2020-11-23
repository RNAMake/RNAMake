#include <rnamake2d/rna_element.h>
#include <iostream>

namespace rnamake2d {
    int RNAElement::num = 0;
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

        auto new_element = RNAElement();
        new_element.type_ = RNAELEMENT::LOOP;

        auto ii = start_index;
        // looping through the specified indices
        while(ii <= end_index) {

            if (pairs_array[ii] < 0) {
                // unpaired region
                new_element.indices_.push_back(ii);
                // checking if the stack is the outermost
                if (last_pair_start >= 0) {
                    auto stack_element = RNAElement();
                    stack_element.type_ = RNAELEMENT::STACK;
                    //stack_element.parent_ = last_parent;
                    if (last_parent != nullptr) {
                        stack_element.parent_id = last_parent->id;
                        //last_parent->children_.push_back(stack_element); // possible problem area
                        last_parent->child_ids.push_back(stack_element.id);
                    }

                    for(auto jj = last_pair_start; jj < ii; ++jj) {
                        stack_element.indices_.push_back(jj);
                        stack_element.indices_.push_back(pairs_array[jj]);
                    }
                    new_element.parent_id = stack_element.id;
                    //new_element.parent_ = &stack_element;
                    stack_element.child_ids.push_back(new_element.id);
                    elements.push_back(std::move(stack_element));
                    last_pair_start = -999;
                    last_pair_end = -999;

                    //stack_element.children_.push_back(new_element);
                }
                ++ii;
            }else if( ii < pairs_array[ii] ) {
                // stack outer pair
                if( ii == start_index and pairs_array[ii] == end_index) {
                    if (last_pair_start < 0 ) {
                        last_pair_start = ii;
                        last_pair_end = pairs_array[ii];
                        last_parent = &new_element;
                    }
                    // reset the parent and then do a recursive call
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
                        //stack_element.parent_ = last_parent;
                        if (last_parent != nullptr ) {
                            stack_element.parent_id = last_parent->id;
                            last_parent->child_ids.push_back(stack_element.id);
                            //last_parent->children_.push_back(stack_element);
                        }

                        for(auto jj = last_pair_start; jj < ii; ++jj) {
                            stack_element.indices_.push_back(jj);
                            stack_element.indices_.push_back(pairs_array[jj]);
                        }
                        new_element.parent_id = stack_element.id;
                        stack_element.child_ids.push_back(new_element.id);
                        elements.push_back(std::move(stack_element));
                        //new_element.parent_ = &stack_element;
                        //stack_element.children_.push_back(new_element);
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
        elements.push_back(std::move(new_element));
    }
    void
    get_rna_elements_from_secstruct_recursive(Ints const& pairs_array,
                                              int start_index,
                                              int end_index,
                                              RNAElems & elements,
                                              int last_pair_start,
                                              int last_pair_end,
                                              std::shared_ptr<RNAElement> last_parent) {

        auto new_element = std::make_shared<RNAElement>();
        new_element->type_ = RNAELEMENT::LOOP;
        /*

        auto ii = start_index;
        // looping through the specified indices
        while(ii <= end_index) {

            if (pairs_array[ii] < 0) {
                // paired region
                new_element->indices_.push_back(ii);
                // checking if the stack is the outermost
                if (last_pair_start >= 0) {
                    auto stack_element = std::make_shared<RNAElement>();
                    stack_element->type_ = RNAELEMENT::STACK;
                    stack_element->parent_ = std::make_shared<RNAElement>(*last_parent);
                    if (last_parent != nullptr) {
                        last_parent->children_.push_back(*stack_element); // possible problem area
                    }

                    for(auto jj = last_pair_start; jj < ii; ++jj) {
                        stack_element->indices_.push_back(jj);
                        stack_element->indices_.push_back(pairs_array[jj]);
                    }

                    elements.push_back(*stack_element);
                    last_pair_start = -999;
                    last_pair_end = -999;
                    new_element->parent_ = std::make_shared<RNAElement>(*last_parent);
                    stack_element->children_.push_back(*new_element);
                }
                ++ii;
            }else if( ii < pairs_array[ii] ) {
                // stack outer pair
                if( ii == start_index and pairs_array[ii] == end_index) {
                    if (last_pair_start < 0 ) {
                        last_pair_start = ii;
                        last_pair_end = pairs_array[ii];
                        last_parent = new_element;
                    }
                    // reset the parent and then do a recursive call
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
                        auto stack_element = std::make_shared<RNAElement>();
                        stack_element->type_ = RNAELEMENT::STACK;
                        stack_element->parent_ = std::make_shared<RNAElement>(*last_parent);
                        if (last_parent != nullptr ) {
                            last_parent->children_.push_back(*stack_element);
                        }

                        for(auto jj = last_pair_start; jj < ii; ++jj) {
                            stack_element->indices_.push_back(jj);
                            stack_element->indices_.push_back(pairs_array[jj]);
                        }
                        elements.push_back(*stack_element);
                        new_element->parent_ = std::make_shared<RNAElement>(*elements.rbegin());
                        stack_element->children_.push_back(*new_element);
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
                            new_element
                    );
                }
                ii = pairs_array[ii] + 1;
            } else {
                LOGE<<start_index<<" "<<end_index;
                LOGE<<"error : should not get here "<<ii<<" "<<pairs_array[ii];
                exit(1);
            }
        }*/
        elements.push_back(*new_element);
    }

    void
    RNAElement::update_family(const std::map<RNAElement *, RNAElement *> &convert) {
//
//        if(parent_ != nullptr) {
//            parent_ = convert.at(parent_.get());
//        }
//
//        for(auto& child : children_) {
//            if(child != nullptr) {
//                child = convert.at(child);
//            }
//        }
    }

    RNAElems
    get_rna_elemnts_from_secstruct(String const& secstruct) {
       auto ptr_elements = std::vector<std::shared_ptr<RNAElement>>{};
       auto pairs_array = get_pairmap_from_secstruct(secstruct);

       std::shared_ptr<RNAElement> dummy(nullptr);

        auto elements = RNAElems();
        get_rna_elements_from_secstruct_recursive(
                pairs_array, 0, secstruct.size() - 1, elements, -999, -999, dummy
                );

        for(auto first = 0; first < elements.size(); ++first )  {
            for(auto second = 0; second < elements.size(); ++second )  {
                if(first == second) {
                    continue;
                }
                auto& elem1 = elements[first];
                auto& elem2 = elements[second];
                if(elem1.parent_id == elem2.id) {
                    if(elem1.parent_ != nullptr ) {
                        std::cout<<"EREERERER"<<std::endl;
                    }
                    elem1.parent_ = &elem2;
                } else if (std::find(elem1.child_ids.begin(), elem1.child_ids.end(), elem2.id) != elem1.child_ids.end()) {
                    auto not_added(true);
                    for(const auto& children : elem1.children_) {
                        if(children->id == elem2.id) {
                            std::cout<<"EREREE"<<std::endl;
                            not_added = false; break;
                        }
                    }
                    if(not_added) {
                        elem1.children_.push_back(&elem2);
                    }
                }
            }
        }

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