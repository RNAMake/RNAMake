//
// Created by Hassan Abdelsamad on 4/30/21.
//

#ifndef RNAMAKE_TESTING_TOOLS_H
#define RNAMAKE_TESTING_TOOLS_H

#include <sstream>

bool compareFile(FILE *fPtr1, FILE *fPtr2) {
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

    // If characters are not same then return -1
    if (ch1 != ch2)
      return -1;

    col += 1;

  } while (ch1 != EOF && ch2 != EOF);

  // If both files have reached end
  if (ch1 == EOF && ch2 == EOF)
    return true;
  else
    return false;
}

#endif // RNAMAKE_TESTING_TOOLS_H
