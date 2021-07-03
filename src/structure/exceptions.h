//
// Created by Joseph Yesselman on 2019-04-24.
//

#ifndef RNAMAKE_NEW_STRUCTURE_EXCEPTIONS_H
#define RNAMAKE_NEW_STRUCTURE_EXCEPTIONS_H

#include "base/types.h"
#include <stdexcept>

namespace structure {

///**
// * Exception for Residues
// */
//class ResidueException : public std::runtime_error {
//public:
//    /**
//     * Standard constructor for ResidueException
//     * @param   message   Error message for residue
//     */
//    ResidueException(String const & message) :
//            std::runtime_error(message) {}
//};

/**
 * Exception for Chains
 */
class ChainException : public std::runtime_error {
public:
    /**
     * Standard constructor for ChainException
     * @param   message   Error message for chain
     */
    ChainException(String const & message) :
            std::runtime_error(message) {}
};

/**
 * Exception for RNA Structure
 */
class RNAStructureException : public std::runtime_error {
public:
    /**
     * Standard constructor for RNAStructureException
     * @param   message   Error message for rna structure
     */
    RNAStructureException(String const & message) :
            std::runtime_error(message) {}
};

}

#endif //RNAMAKE_NEW_STRUCTURE_EXCEPTIONS_H
