//
// Created by Joe Yesselman on 7/7/22.
//

#ifndef RNAMAKE_SRC_RESOURCE_MANAGEMENT_EXCEPTION_HPP_
#define RNAMAKE_SRC_RESOURCE_MANAGEMENT_EXCEPTION_HPP_

#include <base/exception.hpp>

namespace resource_management {

class ResourceManagementException : public std::runtime_error {
public:
  /**
   * Standard constructor for SqliteLibraryException
   * @param   message   Error message for SqliteLibrary
   */
  explicit ResourceManagementException(String const &message)
      : std::runtime_error(message) {}
};

} // namespace resource_management

#endif // RNAMAKE_SRC_RESOURCE_MANAGEMENT_EXCEPTION_HPP_
