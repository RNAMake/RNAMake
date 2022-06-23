//
// Created by Joseph Yesselman on 2019-04-12.
//

#ifndef RNAMAKE_NEW_SEQUENCE_CONSTRAINT_H
#define RNAMAKE_NEW_SEQUENCE_CONSTRAINT_H

#include <secondary_structure/pose.h>
#include <secondary_structure/sequence_tools.h>

namespace secondary_structure {

class SequenceConstraint {
public:
  SequenceConstraint() {}

  virtual ~SequenceConstraint() {}

  virtual SequenceConstraint *clone() const = 0;

public:
  virtual bool violates_constraint(PoseOP p) {
    if (violations(p) > 0) {
      return true;
    } else {
      return false;
    }
  }

  virtual int violations(PoseOP) = 0;
};

class DisallowedSequence : public SequenceConstraint {
public:
  DisallowedSequence(String const &disallowed_sequence) : SequenceConstraint() {
    get_res_types_from_sequence(disallowed_sequence,
                                disallowed_res_type_array_);
  }

  SequenceConstraint *clone() const { return new DisallowedSequence(*this); };

public:
  int violations(PoseOP p) {
    return find_res_types_in_pose(p, disallowed_res_type_array_);
  }

private:
  ResTypes disallowed_res_type_array_;
};

class GCHelixStretchLimit : public SequenceConstraint {
public:
  GCHelixStretchLimit(int length) : SequenceConstraint() { length_ = length; }

  SequenceConstraint *clone() const { return new GCHelixStretchLimit(*this); };

public:
  int violations(PoseOP p) { return find_gc_helix_stretches(p, length_); }

private:
  int length_;
};

typedef std::shared_ptr<SequenceConstraint> SequenceConstraintOP;
typedef std::vector<SequenceConstraintOP> SequenceConstraintOPs;

class SequenceConstraints {
public:
  SequenceConstraints() {}

  ~SequenceConstraints() {}

public:
  void add_sequence_constraint(SequenceConstraintOP seq_constraint) {
    violations_.push_back(0);
    seq_constraints_.push_back(SequenceConstraintOP(seq_constraint->clone()));
  }

  void add_disallowed_sequence(String const &seq) {
    add_sequence_constraint(std::make_shared<DisallowedSequence>(seq));
  }

  void add_gc_helix_stretch_limit(int length) {
    add_sequence_constraint(std::make_shared<GCHelixStretchLimit>(length));
  }

public:
  Ints const &violations(PoseOP p) {
    for (int i = 0; i < seq_constraints_.size(); i++) {
      violations_[i] = seq_constraints_[i]->violations(p);
    }

    return violations_;
  }

  size_t num_constraints() { return seq_constraints_.size(); }

private:
  Ints violations_;
  SequenceConstraintOPs seq_constraints_;
};

} // namespace secondary_structure

#endif // RNAMAKE_NEW_SEQUENCE_CONSTRAINT_H

/*
auto disallowed_sequences = Strings{"AAAA", "CCCC", "GGGG", "UUUU"};
for(auto const & seq : disallowed_sequences) {
auto disallowed_types = secondary_structure::ResTypes();
secondary_structure::get_res_types_from_sequence(seq, disallowed_types);
disallowed_res_type_arrays_.push_back(disallowed_types);

*/