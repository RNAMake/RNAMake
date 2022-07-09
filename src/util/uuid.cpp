//
// Created by Joe Yesselman on 7/6/22.
//

#include <util/uuid.h>

namespace util {
Uuid generate_uuid() { return sole::uuid0(); }

Uuid uuid_from_str(const String &str) { return sole::rebuild(str); }
}