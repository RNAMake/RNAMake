////
//// Created by Joseph Yesselman on 11/12/17.
////
//
//#include <base/assertions.h>
//#include <base/settings.h>
//
// namespace base {
//
//  String
//  get_os_name() {
//#ifdef _WIN32 || _WIN64
//      return  String("windows");
//#elif __unix || __unix__
//      return  String("unix");
//#elif __APPLE__ || __MACH__
//      return String("osx");
//#elif __linux__
//      return String("linux");
//#endif
//      throw std::runtime_error("cannot determine operating system");
//  }
//
//  String
//  base_path() {
//
//      char *base_path;
//      base_path = std::getenv("RNAMAKE_NEW");
//      ensures<std::runtime_error>(
//              base_path != NULL,
//              "cannot find environemental path RNAMAKE_NEW");
//      return String(base_path);
//  }
//
//  String
//  resources_path() {
//      auto base_path_str = base_path();
//      return base_path_str + "/resources/";
//  }
//
//  String
//  unittest_resources_path() {
//      auto base_path_str = base_path();
//      return base_path_str + "/unittests/unittest_resources/";
//  }
//
//
//}