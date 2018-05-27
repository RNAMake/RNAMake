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

        system("rm test.csv");
        system("rm test.bin");

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

        auto path = base_dir() + "/rnamake/lib/RNAMake/unittests/unittest_resources/math/test_2.bin";
        std::ifstream in;
        in.open(path, std::ios::binary);

        auto histo = SixDHistogram(in);
        auto histo_2 = SixDHistogram(in);

        REQUIRE(histo.contains(pA));
        REQUIRE(histo_2.contains(pA));

    }

    SECTION("test read tecto bin file") {
        auto path = base_dir() + "/rnamake/lib/RNAMake/unittests/unittest_resources/math/test_tecto.bin";

        std::ifstream in;
        in.open(path, std::ios::binary);

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

    SECTION("test on tecto data") {
        auto path = base_dir() + "/rnamake/lib/RNAMake/unittests/unittest_resources/math/test.out";

        auto in = std::ifstream();
        in.open(path);
        auto line = String();
        int i = -1;

        auto constraints = std::array<Real2, 6>{
                Real2{-3.5,5.25},
                Real2{-3.5,4.75},
                Real2{-4.00,4.75},
                Real2{9.99999,300},
                Real2{45,335},
                Real2{170,200}};

        auto constraints_2 = std::array<Real2, 6>{
                Real2{-3.5,5.0},
                Real2{-3.5,4.5},
                Real2{-4.0,4.5},
                Real2{5,300},
                Real2{45,330},
                Real2{170,195}};

        auto values = Real6();
        auto bb = BoundingBox(Point(-10, -10, -10), Point(10, 10, 10));
        auto bin_widths = Real6{0.25, 0.25, 0.25, 5, 5, 5};
        auto binner = SixDCoordinateBinner(bb, bin_widths);

        auto count = 0;
        auto count_2 = 0;
        auto fail = 0;
        auto fail_2 = 0;
        while ( in.good() ) {
            i++;
            getline(in, line);
            if (i < 3) { continue; }
            if(line.length() < 5) { break; }
            auto spl = split_str_by_delimiter(line, ",");

            for(int j = 0; j < 6; j++) {
                values[j] = std::stod(spl[j]);
                if(j > 2) { values[j] += 180; }
            }

            auto bin_index = binner.bin_index(values);
            auto bin = binner.bin_from_index(bin_index);
            auto values_2 = binner.bin_to_values(bin);

            fail = 0;
            for(int i = 0; i < 6; i++) {
                if(i != 3 && (constraints[i][0] > values[i] || values[i] > constraints[i][1])) { fail = 1; break; }
                if(i == 3 && (constraints[i][0] < values[i] && values[i] < constraints[i][1])) { fail = 1; break; }
            }

            fail_2 = 0;
            for(int i = 0; i < 6; i++) {
                if(i != 3 && (constraints_2[i][0] > values_2[i] || values_2[i] > constraints_2[i][1])) {
                    fail_2 = 1; break;
                }
                if(i == 3 && (constraints_2[i][0] < values_2[i] && values_2[i] < constraints_2[i][1])) {
                    fail_2 = 1; break;
                }
            }

            if(!fail) { count += 1; }
            if(!fail_2) { count_2 += 1;}

            if(fail != fail_2) {
                std::cout << fail << " " << fail_2 << std::endl;
                for(int i = 0; i < 6; i++) {
                    if (i != 3 && (constraints[i][0] > values[i] || values[i] > constraints[i][1])) {
                        std::cout << i << " " << (constraints[i][0] > values[i]) << " " << (values[i] > constraints[i][1]) << std::endl;
                    }
                    if (i == 3 && (constraints[i][0] < values[i] && values[i] < constraints[i][1])) {
                        std::cout << i << std::endl;
                    }
                }
                std::cout << "bined" << std::endl;

                for(int i = 0; i < 6; i++) {
                    if(i != 3 && (constraints_2[i][0] > values_2[i] || values_2[i] > constraints_2[i][1])) {
                        std::cout << i << " " << (constraints_2[i][0] > values_2[i]) << " " << (values_2[i] > constraints_2[i][1]) << std::endl;                    }
                    if(i == 3 && (constraints_2[i][0] < values_2[i] && values_2[i] < constraints_2[i][1])) {
                        std::cout << i << std::endl;
                    }
                }
                for(auto const & v : values) { std::cout << v << " ";}
                std::cout << std::endl;
                for(auto const & v : values_2) { std::cout << v << " ";}
                std::cout << std::endl;
                exit(0);
            }


        }
        REQUIRE(count == count_2);
    }



}






















