//
// Created by Hassan Abdelsamad on 4/30/21.
//


#ifndef RNAMAKE_TOOLS_H

#include "util/csv.h"


int compareFile(FILE *fPtr1, FILE *fPtr2) {
    char ch1, ch2;
    auto line = 1;
    auto col = 0;

    do {
        // Input character from both files
        ch1 = fgetc(fPtr1);
        ch2 = fgetc(fPtr2);

        // Increment line
        if (ch1 == '\n') {
            line += 1;
            col = 0;
        }

        // If characters are not same then return false
        if (ch1 != ch2)
            return -1;
        col += 1;

    } while (ch1 != EOF && ch2 != EOF);

    /* If both files have reached end */
    if (ch1 == EOF && ch2 == EOF)
        return 1;
    else
        return -1;
}

bool compare_csv(std::string f1, std::string f2) {
    io::LineReader lr1("Logs.csv");
    io::LineReader lr2("Logs_orig.csv");
    auto l1 = lr1.next_line();
    auto l2 = lr2.next_line();
    while (l1 || l2) {
        l1 = lr1.next_line();
        l2 = lr2.next_line();
        if (!(std::strcmp(l1, l2) == 0)) {
            return false;
        }
    }
    return true;
}

#define RNAMAKE_TOOLS_H

#endif //RNAMAKE_TOOLS_H
