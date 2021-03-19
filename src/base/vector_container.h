//
// Created by Joseph Yesselman on 10/23/17.
//

#ifndef RNAMAKE_NEW_VECTOR_CONTAINER_H
#define RNAMAKE_NEW_VECTOR_CONTAINER_H

#include <memory>
#include <exception>

#include <base/assertions.h>
#include <base/types.h>
#include <stdexcept>

namespace base {

// wrapper to avoid doing std::shared_ptr<std::vector<T>> since pybind11 wont allow this
template<typename T>
class VectorContainer {
public:
    inline
    VectorContainer(
            std::vector <T> const & vec) :
            vec_(std::move(vec)) {}

    ~VectorContainer(){}

public:
    inline
    size_t
    size() { return vec_.size(); }

public:
    typedef typename std::vector<T>::const_iterator const_iterator;

    const_iterator
    begin() const noexcept { return vec_.begin(); }

    const_iterator
    end() const noexcept   { return vec_.end(); }


public:

    inline
    T const &
    operator []( Index i ) const {
        expects<std::runtime_error>(
                i < vec_.size(),
                "cannot get element: " + std::to_string(i) + " it is out of range");
        return vec_[i];
    }

    inline
    T const &
    at(Index i) const {
        expects<std::runtime_error>(
                i < vec_.size(),
                "cannot get element: " + std::to_string(i) + " it is out of range");
        return vec_[i];
    }

public:
    inline
    std::vector<T> const &
    get_data() { return vec_; }

private:
    std::vector<T> vec_;

};

template<typename T>
using   VectorContainerUP = std::unique_ptr<VectorContainer<T> >;
template<typename T>
using   VectorContainerOP = std::shared_ptr<VectorContainer<T> >;

}

#endif //RNAMAKE_NEW_VECTOR_CONTAINER_H
