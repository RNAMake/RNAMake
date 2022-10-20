//
// Created by Joseph Yesselman on 2/16/18.
//

#include <random>

#include "../common.hpp"

#include "base/paths.hpp"
#include "math/hashing.h"
#include "math/matrix_3x3.hpp"

TEST_CASE("Test hashing of 6D coords") {
  SUBCASE("test binner from rosetta example") {
    auto lower = math::Vector3(12.5, 16.25, 4.25);
    auto upper = math::Vector3(15.5, 20, 8.5);
    auto bb = math::BoundingBox(lower, upper);

    math::Real6 binwidths;
    binwidths[0] = binwidths[1] = binwidths[2] = 0.25;
    binwidths[3] = binwidths[4] = binwidths[5] = 10;

    auto binner = math::SixDCoordinateBinner(bb, binwidths);

    math::Real6 pA;
    pA[0] = 13.6;
    pA[1] = 19.4;
    pA[2] = 5.3;
    pA[3] = 50;
    pA[4] = 127;
    pA[5] = 76;

    auto bin = binner.bin6(pA);
    auto bin_center = binner.bin_center_point(bin);
    auto bin_index = binner.bin_index(pA);

    auto bin_2 = binner.bin_from_index(bin_index);

    for (int i = 0; i < bin.size(); i++) {
      CHECK(bin[i] == bin_2[i]);
    }
  }
  SUBCASE("test binner wrapping") {
    auto lower = math::Vector3(-5.0, -5.0, -5.0);
    auto upper = math::Vector3(5.0, 5.0, 5.0);
    auto bb = math::BoundingBox(lower, upper);

    math::Real6 binwidths;
    binwidths[0] = binwidths[1] = binwidths[2] = 0.10;
    binwidths[3] = binwidths[4] = binwidths[5] = 5;

    auto binner = math::SixDCoordinateBinner(bb, binwidths);

    math::Real6 pA;
    pA[0] = -4.25;
    pA[1] = 3.42;
    pA[2] = -1.3;
    pA[3] = 360;
    pA[4] = 12.2;
    pA[5] = 2;

    auto bin = binner.bin6(pA);
    auto bin_center = binner.bin_center_point(bin);
    auto bin_index = binner.bin_index(pA);
    auto bin_2 = binner.bin_from_index(bin_index);

    for (int i = 0; i < bin.size(); i++) {
      CHECK(bin[i] == bin_2[i]);
    }
  }
  SUBCASE("test histogram") {
    auto lower = math::Vector3(-5.0, -5.0, -5.0);
    auto upper = math::Vector3(5.0, 5.0, 5.0);
    auto bb = math::BoundingBox(lower, upper);

    math::Real6 binwidths;
    binwidths[0] = binwidths[1] = binwidths[2] = 0.10;
    binwidths[3] = binwidths[4] = binwidths[5] = 5;

    auto histo = math::SixDHistogram(bb, binwidths);

    math::Real6 pA;
    pA[0] = -4.25;
    pA[1] = 3.42;
    pA[2] = -1.3;
    pA[3] = 360;
    pA[4] = 12.2;
    pA[5] = 2;

    histo.add(pA);
    histo.to_text_file("test.csv");

    Strings lines = {};
    base::path::get_lines_from_file("test.csv", lines);
    auto histo_2 = math::SixDHistogram(lines, math::SixDHistogramStrType::TEXT);
    CHECK(histo_2.contains(pA));

    histo.to_binary_file("test.bin");

    std::ifstream in;
    in.open("test.bin", std::ios::binary);
    auto histo_3 = math::SixDHistogram(in);

    CHECK(histo_3.contains(pA));

    system("rm test.csv");
    system("rm test.bin");

    // double test;
    // in.read(reinterpret_cast<char *>(&test), sizeof(test));
    // std::cout << test << std::endl;

    // auto lines2 =base::get_lines_from_file("test.bin");
    // auto histo_3 = math::SixDHistogram(lines2,
    // math::SixDHistogramStrType::BINARY);
  }
  SUBCASE("test two histograms in one file") {
    math::Real6 pA;
    pA[0] = -4.25;
    pA[1] = 3.42;
    pA[2] = -1.3;
    pA[3] = 360;
    pA[4] = 12.2;
    pA[5] = 2;

    auto path = base::path::unittest_resource_path() + "math/test_2.bin";
    std::ifstream in;
    in.open(path, std::ios::binary);

    auto histo = math::SixDHistogram(in);
    auto histo_2 = math::SixDHistogram(in);

    CHECK(histo.contains(pA));
    CHECK(histo_2.contains(pA));
  }
  SUBCASE("test read tecto bin file") {
    auto path = base::path::unittest_resource_path() + "math/test_tecto.bin";

    std::ifstream in;
    in.open(path, std::ios::binary);

    auto histo = math::SixDHistogram(in);
  }
  /*SUBCASE("test histo on tecto data") {
      auto in = std::ifstream();
      in.open("test.out");
      auto line = String();
      int i = -1;

      auto bb = math::BoundingBox(math::Point(-10, -10, -10), math::Point(10,
  10, 10)); auto bin_widths = math::Real6{0.25, 0.25, 0.25, 5, 5, 5}; auto histo
  = math::SixDHistogram(bb, bin_widths); auto values = math::Real6();

      while ( in.good() ) {
          i++;
          getline(in, line);
          if(i == 0) { continue; }
          if(line.length() < 5) { break; }
          auto spl = base::split_str_by_delimiter(line, ",");
          for(int j = 0; j < 6; j++) {
              values[j] = std::stod(spl[j]);
              if(j > 2) { values[j] += 180; }
          }
          histo.add(values);
      }

      histo.to_text_file("test_histo.csv");
  }*/
  SUBCASE("test on tecto data") {
    auto path = base::path::unittest_resource_path() + "math/test.out";

    auto in = std::ifstream();
    in.open(path);
    auto line = String();
    auto i = -1;

    auto constraints = std::array<math::Real2, 6>{
        math::Real2{-3.5, 5.25},  math::Real2{-3.5, 4.75},
        math::Real2{-4.00, 4.75}, math::Real2{9.99999, 300},
        math::Real2{45, 335},     math::Real2{170, 200}};

    auto constraints_2 = std::array<math::Real2, 6>{
        math::Real2{-3.5, 5.0}, math::Real2{-3.5, 4.5}, math::Real2{-4.0, 4.5},
        math::Real2{5, 300},    math::Real2{45, 330},   math::Real2{170, 195}};

    auto values = math::Real6();
    auto bb = math::BoundingBox(math::Vector3(-10, -10, -10),
                                math::Vector3(10, 10, 10));
    auto bin_widths = math::Real6{0.25, 0.25, 0.25, 5, 5, 5};
    auto binner = math::SixDCoordinateBinner(bb, bin_widths);

    auto count = 0;
    auto count_2 = 0;
    auto fail = 0;
    auto fail_2 = 0;
    while (in.good()) {
      i++;
      getline(in, line);
      if (i < 3) {
        continue;
      }
      if (line.length() < 5) {
        break;
      }
      auto spl = base::string::split(line, ",");

      for (int j = 0; j < 6; j++) {
        values[j] = std::stod(spl[j]);
        if (j > 2) {
          values[j] += 180;
        }
      }

      auto bin_index = binner.bin_index(values);
      auto bin = binner.bin_from_index(bin_index);
      auto values_2 = binner.bin_to_values(bin);

      fail = 0;
      for (int i = 0; i < 6; i++) {
        if (i != 3 &&
            (constraints[i][0] > values[i] || values[i] > constraints[i][1])) {
          fail = 1;
          break;
        }
        if (i == 3 &&
            (constraints[i][0] < values[i] && values[i] < constraints[i][1])) {
          fail = 1;
          break;
        }
      }

      fail_2 = 0;
      for (int i = 0; i < 6; i++) {
        if (i != 3 && (constraints_2[i][0] > values_2[i] ||
                       values_2[i] > constraints_2[i][1])) {
          fail_2 = 1;
          break;
        }
        if (i == 3 && (constraints_2[i][0] < values_2[i] &&
                       values_2[i] < constraints_2[i][1])) {
          fail_2 = 1;
          break;
        }
      }

      if (!fail) {
        count += 1;
      }
      if (!fail_2) {
        count_2 += 1;
      }

      if (fail != fail_2) {
        std::cout << fail << " " << fail_2 << std::endl;
        for (int i = 0; i < 6; i++) {
          if (i != 3 && (constraints[i][0] > values[i] ||
                         values[i] > constraints[i][1])) {
            std::cout << i << " " << (constraints[i][0] > values[i]) << " "
                      << (values[i] > constraints[i][1]) << std::endl;
          }
          if (i == 3 && (constraints[i][0] < values[i] &&
                         values[i] < constraints[i][1])) {
            std::cout << i << std::endl;
          }
        }
        std::cout << "bined" << std::endl;

        for (int i = 0; i < 6; i++) {
          if (i != 3 && (constraints_2[i][0] > values_2[i] ||
                         values_2[i] > constraints_2[i][1])) {
            std::cout << i << " " << (constraints_2[i][0] > values_2[i]) << " "
                      << (values_2[i] > constraints_2[i][1]) << std::endl;
          }
          if (i == 3 && (constraints_2[i][0] < values_2[i] &&
                         values_2[i] < constraints_2[i][1])) {
            std::cout << i << std::endl;
          }
        }
        for (auto const &v : values) {
          std::cout << v << " ";
        }
        std::cout << std::endl;
        for (auto const &v : values_2) {
          std::cout << v << " ";
        }
        std::cout << std::endl;
      }
    }
    CHECK(count == count_2);
  }
  SUBCASE("test 3d binner from rosetta example") {
    auto lower = math::Vector3(12.5, 16.25, 4.25);
    auto upper = math::Vector3(15.5, 20, 8.5);
    auto bb = math::BoundingBox(lower, upper);

    math::Real3 binwidths;
    binwidths[0] = binwidths[1] = binwidths[2] = 0.25;

    auto binner = math::ThreeDCoordinateBinner(bb, binwidths);

    auto pA = math::Vector3(13.6, 19.4, 5.3);

    auto bin = binner.bin3(pA);
    auto bin_center = binner.bin_center_point(bin);
    auto bin_index = binner.bin_index(pA);

    auto bin_2 = binner.bin_from_index(bin_index);

    for (int i = 0; i < bin.size(); i++) {
      CHECK(bin[i] == bin_2[i]);
    }
  }
  SUBCASE("test 3d histogram") {
    auto bb = math::BoundingBox(math::Vector3(-100, -100, -100),
                                math::Vector3(100, 100, 100));
    auto bin_widths = math::Real3{0.5, 0.5, 0.5};
    auto histo = math::ThreeDHistogram(bb, bin_widths);

    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    auto points = math::Vector3s();
    for (int i = 0; i < 100; i++) {
      auto p = math::Vector3(dis(gen), dis(gen), dis(gen));
      histo.add(p);
      points.push_back(p);
    }

    for (auto const &p : points) {
      CHECK(histo.contains(p));
    }
  }
  SUBCASE("test BoundingBox with zero width") {
    // should throw an error
    REQUIRE_THROWS_AS(math::BoundingBox(math::Vector3(1, 1, 1), math::Vector3(1, 1, 1)), base::MathException);
  }



}
