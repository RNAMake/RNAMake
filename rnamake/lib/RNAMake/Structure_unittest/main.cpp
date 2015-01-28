//
//  main.cpp
//  Structure_unittest
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "structure.h"
#include "chain.h"

Chain
get_test_chain() {
    String file = "test_str_to_chain.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    input.close();
    Chain c = str_to_chain(line, rts);
    return c;
}

Structure
get_test_structure() {
    String file = "test_str_to_structure.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Structure s = str_to_structure(line, rts);
    return s;
}

int
test_str_to_chain() {
    String file = "test_str_to_chain.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        Chain c = str_to_chain(line, rts);
    }
    return 1;
}

int
test_to_str() {
    String file = "test_str_to_chain.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    input.close();
    Chain c = str_to_chain(line, rts);
    String s = c.to_str();
    Chain c2 = str_to_chain(s, rts);
    //std::cout << c2.residues().size() << std::endl;
    return 1;
}

int
test_subchain() {
    Chain c = get_test_chain();
    Chain sc = c.subchain(1, 5);
    c.to_pdb("original.pdb");
    sc.to_pdb("subchain.pdb");
    
    Residue r1 = c.residues()[1];
    Residue r2 = c.residues()[5];
    Chain sc2 = c.subchain(r1, r2);
    sc2.to_pdb("subchain2.pdb");

    return 1;
}

int
test_to_pdb() {
    Chain c = get_test_chain();
    c.to_pdb("test_to_pdb.pdb");
    return 1;
}

int
test_str_to_structure() {
    String file = "test_str_to_structure.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Structure s = str_to_structure(line, rts);
    return 1;
}

int
test_build_chains() {
    Structure s = get_test_structure();
    Structure s2 = get_test_structure();
    Residues res = s.residues();
    s._build_chains(res);
    for(int i = 0; i < s.residues().size(); i++ ){
        if(s.residues()[i].num() != s2.residues()[i].num()) { return 0; }
    }
    
    return 1;
}

int
test_move() {
    Structure s = get_test_structure();
    Point p (10.0, 0 ,0);
    s.move(p);
    s.chains()[0].to_pdb("moved.pdb");
    return 1;
}

int
test_transform() {
    Structure s = get_test_structure();
    String file = "test_transform.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Matrix r = matrix_from_str(line);
    getline(input, line);
    Point trans = vector_from_str(line);
    getline(input, line);
    Structure s2 = str_to_structure(line, rts);
    Transform t(r, trans);
    s.transform(t);
    return 1;
}


int main(int argc, const char * argv[]) {
    if (test_str_to_chain() == 0)       { std::cout << "test_str_to_chain failed" << std::endl; }
    if (test_to_str() == 0)             { std::cout << "test_to_str failed" << std::endl; }
    if (test_subchain() == 0)           { std::cout << "test_subchain failed" << std::endl; }
    if (test_to_pdb() == 0)             { std::cout << "test_to_pdb failed" << std::endl; }
    if (test_str_to_structure() == 0)   { std::cout << "test_str_to_structure failed" << std::endl; }
    if (test_build_chains() == 0)       { std::cout << "test_build_chains failed" << std::endl; }
    if (test_move() == 0)               { std::cout << "test_move failed" << std::endl; }
    if (test_transform() == 0)          { std::cout << "test_transform failed" << std::endl; }

    return 0;
}











