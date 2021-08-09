//
// Created by Joseph Yesselman on 12/30/17.
//

#include <base/log.h>
#include <util/sqlite/connection.h>


namespace util {
namespace sqlite {

void
create_table(
        Connection & conn,
        TableDetails const & td) {

    auto table_str = "CREATE TABLE " + td.name() + "(";
    auto i = 0;
    for(auto const & col : td) {
        table_str += col.name + " " + col.type;
        i++;
        if(i != td.size()) { table_str += ", "; }
    }
    if(td.has_primary_key()) {
        table_str += ", PRIMARY KEY( ";
        for (auto const & col: td) {
            if(col.is_primary) { table_str += col.name + " "; }
        }
        table_str += ") ";
    }
    table_str += ")";
    conn.exec(table_str);

    LOGI << "creating table: " + table_str + " in database: " + conn.get_database_name();
}

void
insert_many(
        Connection & conn,
        String const & table_name,
        std::vector<Strings> const & data) {

    auto td = conn.get_table_details(table_name);
    auto insert_str = "INSERT INTO " + table_name + "( ";
    auto i = 0;
    for(auto const & col : *td) {
        insert_str += col.name;
        i++;
        if(i != td->size()) { insert_str += ", "; }
    }
    insert_str += ") VALUES ";

    conn.start_transaction();

    auto insert_statement = String("");
    for(auto const & row : data) {
        i = -1;
        insert_statement = insert_str + "(";
        for(auto const & element : row) {
            i++;
            insert_statement += "'" + element + "'";
            if(i != row.size()-1) {
                insert_statement += ",";
            }
        }
        insert_statement += ");";
        conn.exec(insert_statement);
    }

    conn.commit_transaction();
}

}
}