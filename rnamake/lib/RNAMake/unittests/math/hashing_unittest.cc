//
// Created by Joseph Yesselman on 2/16/18.
//


#include "../common.hpp"

#include "base/settings.h"
#include "base/file_io.h"
#include "math/xyz_matrix.h"
#include "math/hashing.h"

TEST_CASE( "Test hashing of 6D coords", "[Hashing]" ) {

    SECTION("test binner from rosetta example") {
        auto lower = Point(12.5, 16.25, 4.25);
        auto upper = Point(15.5, 20, 8.5);
        auto bb = BoundingBox(lower, upper);

        Real6 binwidths;
        binwidths[0] = binwidths[1] = binwidths[2] = 0.25;
        binwidths[3] = binwidths[4] = binwidths[5] = 10;

        auto binner = SixDCoordinateBinner(bb, binwidths);

        Real6 pA;
        pA[0] = 13.6;   pA[1] = 19.4;   pA[2] = 5.3;
        pA[3] = 50;     pA[4] = 127;    pA[5] = 76;

        auto bin = binner.bin6(pA);
        auto bin_center = binner.bin_center_point(bin);
        auto bin_index = binner.bin_index(pA);

        auto bin_2 = binner.bin_from_index(bin_index);

        for (int i = 0; i < bin.size(); i++) {
            REQUIRE(bin[i] == bin_2[i]);
        }
    }

    SECTION("test binner wrapping") {
        auto lower = Point(-5.0, -5.0, -5.0);
        auto upper = Point(5.0, 5.0, 5.0);
        auto bb = BoundingBox(lower, upper);

        Real6 binwidths;
        binwidths[0] = binwidths[1] = binwidths[2] = 0.10;
        binwidths[3] = binwidths[4] = binwidths[5] = 5;

        auto binner = SixDCoordinateBinner(bb, binwidths);

        Real6 pA;
        pA[0] = -4.25;  pA[1] = 3.42;  pA[2] = -1.3;
        pA[3] = 360;    pA[4] = 12.2;  pA[5] = 2;

        auto bin = binner.bin6(pA);
        auto bin_center = binner.bin_center_point(bin);
        auto bin_index = binner.bin_index(pA);
        auto bin_2 = binner.bin_from_index(bin_index);

        for (int i = 0; i < bin.size(); i++) {
            REQUIRE(bin[i] == bin_2[i]);
        }

    }

    SECTION("test histogram") {

        auto lower = Point(-5.0, -5.0, -5.0);
        auto upper = Point(5.0, 5.0, 5.0);
        auto bb = BoundingBox(lower, upper);

        Real6 binwidths;
        binwidths[0] = binwidths[1] = binwidths[2] = 0.10;
        binwidths[3] = binwidths[4] = binwidths[5] = 5;

        auto histo = SixDHistogram(bb, binwidths);

        Real6 pA;
        pA[0] = -4.25;  pA[1] = 3.42;  pA[2] = -1.3;
        pA[3] = 360;    pA[4] = 12.2;  pA[5] = 2;

        histo.add(pA);
        histo.to_text_file("test.csv");

        auto lines = get_lines_from_file("test.csv");
        auto histo_2 = SixDHistogram(lines, SixDHistogramStrType::TEXT);
        REQUIRE(histo_2.contains(pA));

        histo.to_binary_file("test.bin");

        std::ifstream in;
        in.open("test.bin", std::ios::binary);
        auto histo_3 = SixDHistogram(in);

        REQUIRE(histo_3.contains(pA));

        //double test;
        //in.read(reinterpret_cast<char *>(&test), sizeof(test));
        //std::cout << test << std::endl;

        //auto lines2 = get_lines_from_file("test.bin");
        //auto histo_3 = SixDHistogram(lines2, SixDHistogramStrType::BINARY);

    }

    SECTION("test two histograms in one file") {
        Real6 pA;
        pA[0] = -4.25;  pA[1] = 3.42;  pA[2] = -1.3;
        pA[3] = 360;    pA[4] = 12.2;  pA[5] = 2;

        std::ifstream in;
        in.open("test_2.bin", std::ios::binary);

        auto histo = SixDHistogram(in);
        auto histo_2 = SixDHistogram(in);

        REQUIRE(histo.contains(pA));
        REQUIRE(histo_2.contains(pA));

    }

    SECTION("test read tecto bin file") {
        std::ifstream in;
        in.open("test_tecto.bin", std::ios::binary);

        auto histo = SixDHistogram(in);
    }

    /*SECTION("test histo on tecto data") {
        auto in = std::ifstream();
        in.open("test.out");
        auto line = String();
        int i = -1;

        auto bb = BoundingBox(Point(-10, -10, -10), Point(10, 10, 10));
        auto bin_widths = Real6{0.25, 0.25, 0.25, 5, 5, 5};
        auto histo = SixDHistogram(bb, bin_widths);
        auto values = Real6();

        while ( in.good() ) {
            i++;
            getline(in, line);
            if(i == 0) { continue; }
            if(line.length() < 5) { break; }
            auto spl = split_str_by_delimiter(line, ",");
            for(int j = 0; j < 6; j++) {
                values[j] = std::stod(spl[j]);
                if(j > 2) { values[j] += 180; }
            }
            histo.add(values);
        }

        histo.to_text_file("test_histo.csv");
    }*/



}






















