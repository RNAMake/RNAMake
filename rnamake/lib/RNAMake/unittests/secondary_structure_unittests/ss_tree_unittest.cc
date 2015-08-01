//
//  ss_tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "ss_tree_unittest.h"
#include "secondary_structure/ss_tree.h"
#include "secondary_structure/ss_tree_node.h"

int
check_expected(
    sstruct::SS_Tree const & ss_tree,
    Strings const & expected_type,
    Strings const & expected_seq) {
    
    for(auto const & n : ss_tree) {
        if(n->data()->what() != expected_type[n->index()]) {
            return 0;
        }
        if(n->data()->sequence() != expected_seq[n->index()]) {
            return 0;
        }
    }
    return 1;
}


int
SS_TreeUnittest::test_creation() {
    sstruct::SS_Tree ss_tree("GACGUC", "((()))");
    auto expected_type = split_str_by_delimiter("SS_SEQ_BREAK,SS_BP,SS_BP,SS_BP,SS_SEQ_BREAK", ",");
    auto expected_seq  = split_str_by_delimiter("&,G&C,A&U,C&G,&&&", ",");
   
    if(check_expected(ss_tree, expected_type, expected_seq) == 0) {
        return 0;
    }
  
    sstruct::SS_Tree ss_tree1("GCCUUCGGGC", "(((....)))");
    auto expected_type1 = split_str_by_delimiter("SS_SEQ_BREAK,SS_BP,SS_BP,SS_BP,SS_HAIRPIN", ",");
    auto expected_seq1  = split_str_by_delimiter("&,G&C,C&G,C&G,UUCG&", ",");
    
    if(check_expected(ss_tree1, expected_type1, expected_seq1) == 0) {
        return 0;
    }
    
    sstruct::SS_Tree ss_tree2("GGGG+CCCC", "((((+))))");
    auto expected_type2 = split_str_by_delimiter("SS_SEQ_BREAK,SS_BP,SS_BP,SS_BP,SS_BP,SS_SEQ_BREAK", ",");
    auto expected_seq2  = split_str_by_delimiter("&,G&C,G&C,G&C,G&C,&&&", ",");
    
    
    if(check_expected(ss_tree2, expected_type2, expected_seq2) == 0) {
        return 0;
    }
    
    
    sstruct::SS_Tree ss_tree3("GUUGG+CCC", "(..((+)))");
    auto expected_type3 = split_str_by_delimiter("SS_SEQ_BREAK,SS_BP,SS_BULGE,SS_BP,SS_BP,SS_SEQ_BREAK", ",");
    auto expected_seq3  = split_str_by_delimiter("&,G&C,UU&,G&C,G&C,&&&", ",");
    
    if(check_expected(ss_tree3, expected_type3, expected_seq3) == 0) {
        return 0;
    }
    
    sstruct::SS_Tree ss_tree4("GAG&CAG&CAC", "(.(&).(&).)");

    for(auto const & n : ss_tree4) {
     //   std::cout << n->data()->what() << " " << n->data()->sequence() << " " << n->index() << std::endl;
    }
    
    
  
    
    /*String seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA";
    String ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...)))))..";
    */
    //SS_Tree ss_tree3(ss, seq);
    //std::cout << ss_tree3.size() << std::endl;

    //SS_Tree ss_tree4("..(.)..", "AACGGAA");
    
    //SS_Tree ss_tree5("((...(.).(.)...))", "GGAAAGACAGACAAACC");
    //for(auto const & n : ss_tree5) {
    //    std::cout << n.index() << " " << n.data()->what() << " " << n.data()->seq() << std::endl;
    //}

    
    //std::cout << ss_tree.get_node(0)->data()->yb() << std::endl;
    
    return 1;
}


int
SS_TreeUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    return 1;
}