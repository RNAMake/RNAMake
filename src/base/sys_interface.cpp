#include<base/sys_interface.h>



namespace base {

String 
execute_command( const char* cmd ) {
    
    std::array<char, 128> buffer;
    String result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    
    return result;
    
}

String
execute_command( String const & cmd) {
    return execute_command(cmd.c_str());
}

nlohmann::json
execute_command_json( const char* cmd ) {
    return nlohmann::json::parse(execute_command(cmd));
}

nlohmann::json
execute_command_json( String const& cmd )  {
    return nlohmann::json::parse(execute_command(cmd));
}


}
