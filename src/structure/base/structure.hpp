//
// Created by Joe Yesselman on 6/30/22.
//

#ifndef RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_
#define RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_

#include <structure/base/chain.hpp>

namespace structure::base {

template <typename Chain, typename Residue> class Structure {
public:
  typedef std::vector<Residue> Residues;
  typedef std::vector<Chain> Chains;
  typedef std::vector<Chain> ChainsOP;

public:
  /// @brief - constructor
  inline Structure() : _residues(Residues()), _cut_points(Cutpoints()) {}

  inline Structure(Residues &res, Cutpoints &cut_points)
      : _residues(std::move(res)), _cut_points(std::move(cut_points)) {}

  /// @brief - deconstructor
  ~Structure() = default;

public: // res iterator
  typedef typename Residues::const_iterator const_iterator;

  /// @brief - gets the start of the residue chain
  const_iterator begin() const { return _residues.begin(); }

  /// @brief - gets the end of the residue chain
  const_iterator end() const { return _residues.end(); }

public: // get_residue interface
  // TODO Need to remove this function and only use char
  Residue const &get_residue(int num, const String &chain_id,
                             char i_code) const {
    for (auto const &r : _residues) {
      if (num == r.get_num() && chain_id == r.get_chain_id() &&
          i_code == r.get_i_code()) {
        return r;
      }
    }
    throw StructureException("cannot find residue!");

    /*auto ss = std::stringstream();
    ss << "cannot find residue with num: " << num << " chain id: " << chain_id
       << " and i_code";
    throw StructureException(ss.str());   */
  }

  /// @brief - gets a residue by the UUID
  Residue const &get_residue(util::Uuid const &uuid) const {
    for (auto const &r : _residues) {
      if (r.get_uuid() == uuid) {
        return r;
      }
    }
    throw StructureException("cannot find residue by uuid");
  }

  /// @brief - gets a residue by the index
  Residue const &get_residue(Index index) const {
    /*expects<StructureException>(
        index < _residues.size(),
        "cannot get residue " + std::to_string(index) + " only " +
            std::to_string(_residues.size()) + " total residues");   */

    return _residues[index];
  }

  /// @brief - gets the index of a residue
  int get_res_index(Residue const &res) const {
    int i = -1;
    for (auto const &r : _residues) {
      i++;
      if (r == res) {
        return i;
      }
    }
    throw StructureException("cannot find index for res: " +
                             std::to_string(res.get_num()));
  }

public: // getters /////////////////////////////////////////////////////////////
  /// @brief -
  ChainsOP get_chains() const {
    auto pos = 0;
    auto res = Residues();
    auto chains = std::vector<Chain>();
    auto i = 0;
    for (auto const &r : _residues) {
      if (_cut_points[pos] == i) {
        auto c = Chain(res);
        chains.push_back(c);
        res = Residues{Residue(r)};
        pos += 1;
      } else {
        res.push_back(Residue(r));
      }
      i++;
    }
    if (res.size() > 0) {
      chains.push_back(Chain(res));
    }
    return chains;
  }

  /// @brief - retuns cutpoints
  [[nodiscard]] const Cutpoints &get_cutpoints() const { return _cut_points; }

  /// @brief - counts and returns the number of residues
  [[nodiscard]] size_t get_num_residues() const { return _residues.size(); }

  /// @brief - returns residues
  [[nodiscard]] const Residues get_residues() const { return _residues; }

  /// @brief - counts and returns the number of chains
  [[nodiscard]] size_t get_num_chains() const { return _cut_points.size(); }

  /// @brief -
  [[nodiscard]] String get_sequence() const {
    auto i = -1;
    auto seq = String("");
    auto pos = 0;
    for (auto const &r : _residues) {
      i++;
      if (_cut_points[pos] == i) {
        seq += "&";
        pos++;
      }
      seq += r.get_name();
    }
    return seq;
  }

  /// @brief - get string representation
  [[nodiscard]] String get_str() const {
    String str;
    for(auto const & c : get_chains()) {
        str += c.get_str() + ":";
    }
    return str;
  }

  /// @brief - checks if the residue is the start of a chain
  bool is_residue_start_of_chain(Residue const &r) const {
    auto res_index = get_res_index(r);
    if (res_index == 0) {
      return true;
    }
    for (auto const c : _cut_points) {
      if (res_index == c) {
        return true;
      }
    }
    return false;
  }

  /// @brief - checks if the residue is the end of a chain
  bool is_residue_end_of_chain(Residue const &r) const {
    auto res_index = get_res_index(r);
    for (auto const c : _cut_points) {
      if (res_index == c - 1) {
        return true;
      }
    }
    return false;
  }

public:
  /// @brief - moves each residue in a chain by the specified position vector
  void move(const math::Vector3 &p) {
    for (auto &r : _residues) {
      r.move(p);
    }
  }

  /// @brief - rotates each residue in a chain by the specified rotation matrix
  void rotate(const math::Matrix3x3 &rot) {
    for (auto &r : _residues) {
      r.rotate(rot);
    }
  }


  inline bool compare_structures(const Structure &structure) const {
    std::vector<std::vector<std::string>> this_atom_coords;
    for(auto const &c : get_chains()) {
      auto chain_str = ::base::string::split(c.get_str(), ",");
      for (auto chain_data : chain_str) {
        auto chain_datum = ::base::string::split(chain_data, " ");
        if (chain_datum.size() > 1) {
          std::vector<std::string> coords;
          coords.push_back(chain_datum[0]);
          coords.push_back(chain_datum[1]);
          coords.push_back(chain_datum[2]);
          coords.push_back(chain_datum[3]);
          this_atom_coords.push_back(coords);
        }
      }
    }

    auto other_chains = structure.get_chains();
    std::vector<std::vector<std::string>> other_atom_coords;
    for (auto const &c : other_chains) {
      auto chain_str = ::base::string::split(c.get_str(), ",");
      for (auto chain_data : chain_str) {
        auto chain_datum = ::base::string::split(chain_data, " ");
        if (chain_datum.size() > 1) {
          std::vector<std::string> coords;
          coords.push_back(chain_datum[0]);
          coords.push_back(chain_datum[1]);
          coords.push_back(chain_datum[2]);
          coords.push_back(chain_datum[3]);
          other_atom_coords.push_back(coords);
        }
      }
    }

    if (this_atom_coords.size() != other_atom_coords.size()) {
      std::cout << "Coordinate size mismatch\n";
      return false;
    }

    const float threshold = 0.1; // Tweak this?
    for (int i = 0; i < this_atom_coords.size(); i++) {
      auto these = this_atom_coords[i];
      auto those = other_atom_coords[i];

      // Check atom name
      if (these[0] != those[0]) {
        std::cout << "Atom mismatch in structure: " << these[0] << " != " << those[0] << "\n";
        return false;
      }

      // Check x coordinates
      float x_diff = std::stof(these[1]) - std::stof(those[1]);
      if (abs(x_diff) > threshold) {
        std::cout << these[0] << " X Ccoordinate mismatch: " << these[1] << " - " << those[1] << " > threshold\n";
      }

      // Check y coordinates
      float y_diff = std::stof(these[2]) - std::stof(those[2]);
      if (abs(y_diff) > threshold) {
        std::cout << these[0] << " Y Ccoordinate mismatch: " << these[2] << " - " << those[2] << " > threshold\n";
      }

      // Check z coordinates
      float z_diff = std::stof(these[3]) - std::stof(those[3]);
      if (abs(z_diff) > threshold) {
        std::cout << these[0] << " Z Ccoordinate mismatch: " << these[3] << " - " << those[3] << " > threshold\n";
      }
    }
    return true;
  }

  inline bool operator==(const Structure &structure) const {
    if (compare_structures(structure)) {
      return true;
    }
    return false;
  }

  inline bool operator!=(const Structure &structure) const {
    if (compare_structures(structure)) {
      return false;
    }
    return true;
  }

private:
  Residues _residues;
  Cutpoints _cut_points;
};

} // namespace structure::base

#endif // RNAMAKE_SRC_STRUCTURE_BASE_STRUCTURE_HPP_
