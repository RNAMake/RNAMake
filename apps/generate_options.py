import csv
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--app', type=str, help='name of the app')


args = parser.parse_args()

app_name = args.app

if not (os.path.isdir(app_name) and os.path.isfile(app_name + '/options.csv')):
        quit('No file')

file = open(args.app + '/' + args.app + '_options.h', 'w+')

s = '''
#include "design_rna_scaffold/design_rna_scaffold.h"

#ifndef ARGS_H
#define ARGS_H

String
valid_pdb (String &path) {
    auto ending = path.substr(path.size() - 4);
    return ending == ".pdb" ? String{""} : String{"the file specified by --pdb must end in .pdb"};
}

String
valid_bp (String &bp) {
    const auto bp_pattern = std::regex("\\b[A-Z][0-9]*-[A-Z][0-9]*\\b");
    auto sm = std::smatch{};
    std::regex_match(bp, sm, bp_pattern);

    return sm.size() == 1 ? String{""} : String{bp + " is an invalid bp format"};
}

struct Parameters {

'''

csv_file = app_name + '/options.csv'
df = pd.read_csv(csv_file, quoting=csv.QUOTE_NONE)

group_list = {}
struct_list = {}

for index, row in df.iterrows():
    if row['group'] not in group_list.keys():
        group_list[row['group']] = f"\tstruct {row['struct']}{{\n"
        struct_list[row['group']] = row['struct']


for index, row in df.iterrows():
    group_list[row['group']] += f"\t\t{row['type']} {row['parameter']} = {row['default']};\n"

for group in group_list:
    group_list[group] += '\t};\n\n'
    s += group_list[group]

for group in group_list:
    s += f"\n\t{struct_list[group]} {group} = {struct_list[group]}();"

s += '''

};

Parameters parameters_ = Parameters();

void
DesignRNAScaffold::setup_options () {

'''


for group in group_list:
    s += f"\tapp_.add_option_group(\"{group}\");\n"

for index, row in df.iterrows():
    if row['flag']:
        s += f"\tapp_.add_flag(\"--{row['name']}\", parameters_.{row['group']}.{row['parameter']}, \"{row['help']}\")"
    else:
        s += f"\tapp_.add_option(\"--{row['name']}\", parameters_.{row['group']}.{row['parameter']}, \"{row['help']}\")"
    if not pd.isnull(row['default']):
        s += f"\n\t\t->default_val({row['default']})"

    s += f"\n\t\t->group(\"{row['group']}\")"

    s += ';\n\n'

s += '''
}

#endif //ARGS_H
'''

file.write(s)
