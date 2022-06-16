//
// Created by Joseph Yesselman on 5/15/18.
//

#ifndef RNAMAKE_NEW_SEGMENT_GRAPH_FWD_H_H
#define RNAMAKE_NEW_SEGMENT_GRAPH_FWD_H_H

#include <memory>

template <typename SegmentType, typename AlignerType> class SegmentGraph;

template <typename SegmentType, typename AlignerType>
using SegmentGraphOP = std::shared_ptr<SegmentGraph<SegmentType, AlignerType>>;

#endif // RNAMAKE_NEW_SEGMENT_GRAPH_FWD_H_H
