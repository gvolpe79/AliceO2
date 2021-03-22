// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// This file has been generated. Do not modify it by hand or your changes might
// be lost.
//
#include "CathodeSegmentationCreator.h"

namespace o2
{
namespace mch
{
namespace mapping
{
namespace impl4
{
CathodeSegmentation* createSegType2(bool isBendingPlane)
{
  if (isBendingPlane) {
    return new CathodeSegmentation{
      2,
      true,
      /* PG */
      {{1, 0, 0, 50, 2}, {2, 12, 0, 45, 4}, {3, 12, 0, 40, 4}, {6, 8, 0, 35, 0}, {7, 15, 0, 30, 0}, {8, 13, 0, 27.5, 4}, {9, 14, 0, 22.5, 0}, {10, 7, 0, 20, 0}, {11, 8, 0, 15, 0}, {12, 15, 0, 10, 0}, {13, 13, 0, 7.5, 4}, {14, 14, 0, 2.5, 0}, {15, 7, 0, 0, 0}, {104, 8, 1, -50, 0}, {105, 15, 1, -60, 0}, {106, 13, 1, -65, 4}, {107, 18, 1, -75, 0}, {111, 8, 1, -10, 0}, {112, 15, 1, -20, 0}, {113, 13, 1, -25, 4}, {114, 14, 1, -35, 0}, {115, 7, 1, -40, 0}, {201, 5, 1, -75, -20}, {202, 6, 1, -70, -20}, {203, 11, 1, -65, -20}, {204, 17, 1, -60, -20}, {205, 10, 1, -50, -20}, {209, 9, 1, -40, -20}, {210, 16, 1, -35, -20}, {211, 11, 1, -25, -20}, {212, 17, 1, -20, -20}, {213, 10, 1, -10, -20}, {304, 1, 0, 40, -20}, {305, 2, 0, 42.5, -20}, {306, 3, 0, 45, -20}, {307, 4, 0, 50, -20}, {315, 9, 0, 0, -20}, {316, 16, 0, 2.5, -20}, {317, 11, 0, 7.5, -20}, {318, 17, 0, 10, -20}, {319, 10, 0, 15, -20}, {320, 9, 0, 20, -20}, {321, 16, 0, 22.5, -20}, {322, 11, 0, 27.5, -20}, {323, 17, 0, 30, -20}, {324, 10, 0, 35, -20}},
      /* PGT */
      {/* C10 */ {3, 36, {28, -1, -1, 29, -1, -1, 30, -1, -1, 31, -1, -1, 58, -1, -1, 55, -1, -1, 54, -1, -1, 52, -1, -1, 49, -1, -1, 48, -1, -1, 46, -1, -1, 43, -1, -1, 42, -1, -1, 40, -1, -1, 39, -1, -1, 32, 4, -1, 37, 0, -1, 34, 3, -1, 33, 1, -1, 36, 6, -1, 35, 2, -1, 38, 5, -1, 41, 8, -1, 44, 7, -1, 45, 13, -1, 47, 11, -1, 50, 14, -1, 51, 12, -1, 53, 19, -1, 56, 16, 25, 57, 21, 22, 59, 20, 18, 60, 23, 17, 61, 24, 15, 62, 26, 10, 63, 27, 9}},
       /* C6 */ {2, 48, {25, 35, 22, 38, 18, 41, 17, 44, 15, 45, 10, 47, 9, 50, 4, 51, 0, 53, 3, 56, 1, 57, 6, 59, 2, 60, 5, 61, 8, 62, 7, 63, 13, -1, 11, -1, 14, -1, 12, -1, 19, -1, 16, -1, 21, -1, 20, -1, 23, -1, 24, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1, 58, -1, 55, -1, 54, -1, 52, -1, 49, -1, 48, -1, 46, -1, 43, -1, 42, -1, 40, -1, 39, -1, 32, -1, 37, -1, 34, -1, 33, -1, 36, -1}},
       /* C7 */ {2, 48, {-1, 58, -1, 55, -1, 54, -1, 52, -1, 49, -1, 48, -1, 46, -1, 43, -1, 42, -1, 40, -1, 39, -1, 32, -1, 37, -1, 34, -1, 33, -1, 36, 25, 35, 22, 38, 18, 41, 17, 44, 15, 45, 10, 47, 9, 50, 4, 51, 0, 53, 3, 56, 1, 57, 6, 59, 2, 60, 5, 61, 8, 62, 7, 63, 13, -1, 11, -1, 14, -1, 12, -1, 19, -1, 16, -1, 21, -1, 20, -1, 23, -1, 24, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1}},
       /* C8 */ {2, 48, {-1, 13, -1, 11, -1, 14, -1, 12, -1, 19, -1, 16, -1, 21, -1, 20, -1, 23, -1, 24, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1, 58, -1, 55, -1, 54, -1, 52, -1, 49, -1, 48, -1, 46, -1, 43, -1, 42, -1, 40, -1, 39, -1, 32, -1, 37, -1, 34, -1, 33, -1, 36, 25, 35, 22, 38, 18, 41, 17, 44, 15, 45, 10, 47, 9, 50, 4, 51, 0, 53, 3, 56, 1, 57, 6, 59, 2, 60, 5, 61, 8, 62, 7, 63}},
       /* C9 */ {3, 36, {25, 49, 56, 22, 48, 57, 18, 46, 59, 17, 43, 60, 15, 42, 61, 10, 40, 62, 9, 39, 63, 4, 32, -1, 0, 37, -1, 3, 34, -1, 1, 33, -1, 6, 36, -1, 2, 35, -1, 5, 38, -1, 8, 41, -1, 7, 44, -1, 13, 45, -1, 11, 47, -1, 14, 50, -1, 12, 51, -1, 19, 53, -1, 16, -1, -1, 21, -1, -1, 20, -1, -1, 23, -1, -1, 24, -1, -1, 26, -1, -1, 27, -1, -1, 28, -1, -1, 29, -1, -1, 30, -1, -1, 31, -1, -1, 58, -1, -1, 55, -1, -1, 54, -1, -1, 52, -1, -1}},
       /* I1 */ {1, 64, {25, 22, 18, 17, 15, 10, 9, 4, 0, 3, 1, 6, 2, 5, 8, 7, 13, 11, 14, 12, 19, 16, 21, 20, 23, 24, 26, 27, 28, 29, 30, 31, 58, 55, 54, 52, 49, 48, 46, 43, 42, 40, 39, 32, 37, 34, 33, 36, 35, 38, 41, 44, 45, 47, 50, 51, 53, 56, 57, 59, 60, 61, 62, 63}},
       /* L18 */ {2, 40, {13, -1, 11, -1, 14, -1, 12, -1, 19, -1, 16, -1, 21, -1, 20, -1, 23, -1, 24, -1, 26, -1, 27, -1, 28, -1, 29, -1, 30, -1, 31, -1, 58, -1, 55, -1, 54, -1, 52, -1, 49, -1, 48, -1, 46, -1, 43, -1, 42, -1, 40, -1, 39, -1, 32, -1, 37, -1, 34, -1, 33, -1, 36, -1, 35, 63, 38, 62, 41, 61, 44, 60, 45, 59, 47, 57, 50, 56, 51, 53}},
       /* L5 */ {2, 40, {23, 20, 24, 21, 26, 16, 27, 19, 28, 12, 29, 14, 30, 11, 31, 13, 58, 7, 55, 8, 54, 5, 52, 2, 49, 6, 48, 1, 46, 3, 43, 0, 42, 4, 40, 9, 39, 10, 32, 15, 37, 17, 34, 18, 33, 22, 36, 25, 35, -1, 38, -1, 41, -1, 44, -1, 45, -1, 47, -1, 50, -1, 51, -1, 53, -1, 56, -1, 57, -1, 59, -1, 60, -1, 61, -1, 62, -1, 63, -1}},
       /* L6 */ {2, 40, {42, 43, 40, 46, 39, 48, 32, 49, 37, 52, 34, 54, 33, 55, 36, 58, 35, 31, 38, 30, 41, 29, 44, 28, 45, 27, 47, 26, 50, 24, 51, 23, 53, 20, 56, 21, 57, 16, 59, 19, 60, 12, 61, 14, 62, 11, 63, 13, -1, 7, -1, 8, -1, 5, -1, 2, -1, 6, -1, 1, -1, 3, -1, 0, -1, 4, -1, 9, -1, 10, -1, 15, -1, 17, -1, 18, -1, 22, -1, 25}},
       /* L7 */ {2, 40, {25, -1, 22, -1, 18, -1, 17, -1, 15, -1, 10, -1, 9, -1, 4, -1, 0, -1, 3, -1, 1, -1, 6, -1, 2, -1, 5, -1, 8, -1, 7, -1, 13, 63, 11, 62, 14, 61, 12, 60, 19, 59, 16, 57, 21, 56, 20, 53, 23, 51, 24, 50, 26, 47, 27, 45, 28, 44, 29, 41, 30, 38, 31, 35, 58, 36, 55, 33, 54, 34, 52, 37, 49, 32, 48, 39, 46, 40, 43, 42}},
       /* L8 */ {2, 40, {-1, 63, -1, 62, -1, 61, -1, 60, -1, 59, -1, 57, -1, 56, -1, 53, -1, 51, -1, 50, -1, 47, -1, 45, -1, 44, -1, 41, -1, 38, -1, 35, 25, 36, 22, 33, 18, 34, 17, 37, 15, 32, 10, 39, 9, 40, 4, 42, 0, 43, 3, 46, 1, 48, 6, 49, 2, 52, 5, 54, 8, 55, 7, 58, 13, 31, 11, 30, 14, 29, 12, 28, 19, 27, 16, 26, 21, 24, 20, 23}},
       /* O10 */ {2, 32, {31, 58, 30, 55, 29, 54, 28, 52, 27, 49, 26, 48, 24, 46, 23, 43, 20, 42, 21, 40, 16, 39, 19, 32, 12, 37, 14, 34, 11, 33, 13, 36, 7, 35, 8, 38, 5, 41, 2, 44, 6, 45, 1, 47, 3, 50, 0, 51, 4, 53, 9, 56, 10, 57, 15, 59, 17, 60, 18, 61, 22, 62, 25, 63}},
       /* O25 */ {2, 32, {58, 25, 55, 22, 54, 18, 52, 17, 49, 15, 48, 10, 46, 9, 43, 4, 42, 0, 40, 3, 39, 1, 32, 6, 37, 2, 34, 5, 33, 8, 36, 7, 35, 13, 38, 11, 41, 14, 44, 12, 45, 19, 47, 16, 50, 21, 51, 20, 53, 23, 56, 24, 57, 26, 59, 27, 60, 28, 61, 29, 62, 30, 63, 31}},
       /* O9 */ {2, 32, {63, 25, 62, 22, 61, 18, 60, 17, 59, 15, 57, 10, 56, 9, 53, 4, 51, 0, 50, 3, 47, 1, 45, 6, 44, 2, 41, 5, 38, 8, 35, 7, 36, 13, 33, 11, 34, 14, 37, 12, 32, 19, 39, 16, 40, 21, 42, 20, 43, 23, 46, 24, 48, 26, 49, 27, 52, 28, 54, 29, 55, 30, 58, 31}},
       /* Z1 */ {3, 40, {-1, 0, 4, -1, 3, 9, -1, 1, 10, -1, 6, 15, -1, 2, 17, -1, 5, 18, -1, 8, 22, -1, 7, 25, -1, 13, -1, -1, 11, -1, -1, 14, -1, -1, 12, -1, -1, 19, -1, -1, 16, -1, -1, 21, -1, -1, 20, -1, -1, 23, -1, -1, 24, -1, -1, 26, -1, -1, 27, -1, -1, 28, -1, -1, 29, -1, -1, 30, -1, -1, 31, -1, 63, 58, -1, 62, 55, -1, 61, 54, -1, 60, 52, -1, 59, 49, -1, 57, 48, -1, 56, 46, -1, 53, 43, -1, 51, 42, -1, 50, 40, -1, 47, 39, -1, 45, 32, -1, 44, 37, -1, 41, 34, -1, 38, 33, -1, 35, 36, -1}},
       /* Z2 */ {3, 40, {53, 51, -1, 56, 50, -1, 57, 47, -1, 59, 45, -1, 60, 44, -1, 61, 41, -1, 62, 38, -1, 63, 35, -1, -1, 36, -1, -1, 33, -1, -1, 34, -1, -1, 37, -1, -1, 32, -1, -1, 39, -1, -1, 40, -1, -1, 42, -1, -1, 43, -1, -1, 46, -1, -1, 48, -1, -1, 49, -1, -1, 52, -1, -1, 54, -1, -1, 55, -1, -1, 58, -1, -1, 31, 25, -1, 30, 22, -1, 29, 18, -1, 28, 17, -1, 27, 15, -1, 26, 10, -1, 24, 9, -1, 23, 4, -1, 20, 0, -1, 21, 3, -1, 16, 1, -1, 19, 6, -1, 12, 2, -1, 14, 5, -1, 11, 8, -1, 13, 7}},
       /* Z3 */ {3, 40, {7, 13, -1, 8, 11, -1, 5, 14, -1, 2, 12, -1, 6, 19, -1, 1, 16, -1, 3, 21, -1, 0, 20, -1, 4, 23, -1, 9, 24, -1, 10, 26, -1, 15, 27, -1, 17, 28, -1, 18, 29, -1, 22, 30, -1, 25, 31, -1, -1, 58, -1, -1, 55, -1, -1, 54, -1, -1, 52, -1, -1, 49, -1, -1, 48, -1, -1, 46, -1, -1, 43, -1, -1, 42, -1, -1, 40, -1, -1, 39, -1, -1, 32, -1, -1, 37, -1, -1, 34, -1, -1, 33, -1, -1, 36, -1, -1, 35, 63, -1, 38, 62, -1, 41, 61, -1, 44, 60, -1, 45, 59, -1, 47, 57, -1, 50, 56, -1, 51, 53}},
       /* Z4 */ {3, 40, {-1, 36, 35, -1, 33, 38, -1, 34, 41, -1, 37, 44, -1, 32, 45, -1, 39, 47, -1, 40, 50, -1, 42, 51, -1, 43, 53, -1, 46, 56, -1, 48, 57, -1, 49, 59, -1, 52, 60, -1, 54, 61, -1, 55, 62, -1, 58, 63, -1, 31, -1, -1, 30, -1, -1, 29, -1, -1, 28, -1, -1, 27, -1, -1, 26, -1, -1, 24, -1, -1, 23, -1, -1, 20, -1, -1, 21, -1, -1, 16, -1, -1, 19, -1, -1, 12, -1, -1, 14, -1, -1, 11, -1, -1, 13, -1, 25, 7, -1, 22, 8, -1, 18, 5, -1, 17, 2, -1, 15, 6, -1, 10, 1, -1, 9, 3, -1, 4, 0, -1}},
       /* Z5 */
       {3,
        40,
        {-1, 0, 4, -1, 3, 9, -1, 1, 10, -1, 6, 15, -1, 2, 17,
         -1, 5, 18, -1, 8, 22, -1, 7, 25, -1, 13, -1, -1, 11, -1,
         -1, 14, -1, -1, 12, -1, -1, 19, -1, -1, 16, -1, -1, 21, -1,
         -1, 20, -1, -1, 23, -1, -1, 24, -1, -1, 26, -1, -1, 27, -1,
         -1, 28, -1, -1, 29, -1, -1, 30, -1, -1, 31, -1, 63, 58, -1,
         62, 55, -1, 61, 54, -1, 60, 52, -1, 59, 49, -1, 57, 48, -1,
         56, 46, -1, 53, 43, -1, 51, 42, -1, 50, 40, -1, 47, 39, -1,
         45, 32, -1, 44, 37, -1, 41, 34, -1, 38, 33, -1, 35, 36, -1}}},
      /* PS */
      {{2.5, 0.5}, {5, 0.5}}};
  } else {
    return new CathodeSegmentation{
      2,
      false,
      /* PG */
      {{1028, 4, 0, 40, 2.5},
       {1029, 3, 0, 47.1428566, 5},
       {1040, 7, 0, -3.996802889e-15, 0},
       {1041, 7, 0, 5.714285851, 0},
       {1042, 7, 0, 11.4285717, 0},
       {1043, 7, 0, 17.1428566, 0},
       {1044, 7, 0, 22.8571434, 0},
       {1045, 7, 0, 28.5714283, 0},
       {1046, 7, 0, 34.2857132, 0},
       {1125, 10, 1, -74.2857132, 0},
       {1126, 10, 1, -62.8571434, 0},
       {1127, 10, 1, -51.42856979, 0},
       {1132, 5, 1, -40, 0},
       {1133, 10, 1, -25.7142849, 0},
       {1134, 6, 1, -14.28571415, 0},
       {1230, 9, 1, -51.42856979, -20},
       {1231, 9, 1, -62.8571434, -20},
       {1232, 9, 1, -74.2857132, -20},
       {1238, 12, 1, -10, -20},
       {1239, 14, 1, -20, -20},
       {1240, 13, 1, -31.4285717, -20},
       {1241, 11, 1, -40, -20},
       {1325, 2, 0, 49.2857132, -20},
       {1326, 1, 0, 44.2857132, -20},
       {1327, 0, 0, 40, -20},
       {1332, 8, 0, 34.2857132, -20},
       {1333, 8, 0, 28.5714283, -20},
       {1334, 8, 0, 22.8571434, -20},
       {1335, 8, 0, 17.1428566, -20},
       {1336, 8, 0, 11.4285717, -20},
       {1337, 8, 0, 5.714285851, -20},
       {1338, 8, 0, -3.996802889e-15, -20}},
      /* PGT */
      {/* C1 */ {7, 10, {51, 33, 49, 26, 13, 9, -1, 53, 36, 48, 27, 11, 4, -1, 56, 35, 46, 28, 14, 0, -1, 57, 38, 43, 29, 12, 3, -1, 59, 41, 42, 30, 19, 1, 25, 60, 44, 40, 31, 16, 6, 22, 61, 45, 39, 58, 21, 2, 18, 62, 47, 32, 55, 20, 5, 17, 63, 50, 37, 54, 23, 8, 15, -1, -1, 34, 52, 24, 7, 10}},
       /* C2 */ {7, 10, {60, 41, 42, 30, 19, 1, 25, 61, 44, 40, 31, 16, 6, 22, 62, 45, 39, 58, 21, 2, 18, 63, 47, 32, 55, 20, 5, 17, -1, 50, 37, 54, 23, 8, 15, -1, 51, 34, 52, 24, 7, 10, -1, 53, 33, 49, 26, 13, 9, -1, 56, 36, 48, 27, 11, 4, -1, 57, 35, 46, 28, 14, 0, -1, 59, 38, 43, 29, 12, 3}},
       /* C3 */
       {13, 10, {50, 37, 54, 23, 14, 8, 1, 4, 10, 17, 18, 22, 25, 51, 34, 52, 24, 12, 7, 6, 0, 9, 15, -1, -1, -1, 53, 33, 49, 26, 19, 13, 2, 3, -1, -1, -1, -1, -1, 56, 36, 48, 27, 16, 11, 5, -1, -1, -1, -1, -1, -1, 57, 35, 46, 28, 21, -1, -1, -1, -1, -1, -1, -1, -1, 59, 38, 43, 29, 20, -1, -1, -1, -1, -1, -1, -1, -1, 60, 41, 42, 30, -1, -1, -1, -1, -1, -1, -1, -1, -1, 61, 44, 40, 31, -1, -1, -1, -1, -1, -1, -1, -1, -1, 62, 45, 39, 58, -1, -1, -1, -1, -1, -1, -1, -1, -1, 63, 47, 32, 55, -1, -1, -1, -1, -1, -1, -1, -1, -1}},
       /* C4 */ {16, 6, {-1, 15, 1, 13, 21, 28, 54, 42, -1, -1, -1, -1, -1, -1, -1, -1, -1, 10, 6, 11, 20, 29, 52, 40, -1, -1, -1, -1, -1, -1, -1, -1, 25, 9, 2, 14, 23, 30, 49, 39, 33, 41, -1, -1, -1, -1, -1, -1, 22, 4, 5, 12, 24, 31, 48, 32, 36, 44, 50, -1, -1, -1, -1, -1, 18, 0, 8, 19, 26, 58, 46, 37, 35, 45, 51, 56, 59, -1, -1, -1, 17, 3, 7, 16, 27, 55, 43, 34, 38, 47, 53, 57, 60, 61, 62, 63}},
       /* C5 */ {11, 7, {25, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, 22, 0, 8, 19, 26, 58, 46, 37, 41, 53, 62, 18, 3, 7, 16, 27, 55, 43, 34, 44, 56, 63, 17, 1, 13, 21, 28, 54, 42, 33, 45, 57, -1, 15, 6, 11, 20, 29, 52, 40, 36, 47, 59, -1, 10, 2, 14, 23, 30, 49, 39, 35, 50, 60, -1, 9, 5, 12, 24, 31, 48, 32, 38, 51, 61, -1}},
       /* L3 */ {20, 4, {17, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 18, 9, 1, 8, 14, 16, 23, 27, 30, 55, 49, 43, 39, 34, 35, 44, 50, 56, 60, 63, 22, 10, 3, 5, 11, 19, 20, 26, 29, 58, 52, 46, 40, 37, 36, 41, 47, 53, 59, 62, 25, 15, 0, 2, 13, 12, 21, 24, 28, 31, 54, 48, 42, 32, 33, 38, 45, 51, 57, 61}},
       /* L4 */ {20, 4, {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 44, 51, 59, 63, 18, 10, 0, 6, 8, 11, 19, 20, 26, 29, 58, 52, 46, 40, 37, 36, 41, 50, 57, 62, 22, 15, 4, 1, 5, 13, 12, 21, 24, 28, 31, 54, 48, 42, 32, 33, 38, 47, 56, 61, 25, 17, 9, 3, 2, 7, 14, 16, 23, 27, 30, 55, 49, 43, 39, 34, 35, 45, 53, 60}},
       /* O1 */ {8, 8, {4, 7, 20, 31, 43, 36, 51, 63, 9, 8, 21, 30, 46, 33, 50, 62, 10, 5, 16, 29, 48, 34, 47, 61, 15, 2, 19, 28, 49, 37, 45, 60, 17, 6, 12, 27, 52, 32, 44, 59, 18, 1, 14, 26, 54, 39, 41, 57, 22, 3, 11, 24, 55, 40, 38, 56, 25, 0, 13, 23, 58, 42, 35, 53}},
       /* O2 */ {8, 8, {53, 35, 42, 58, 23, 13, 0, 25, 56, 38, 40, 55, 24, 11, 3, 22, 57, 41, 39, 54, 26, 14, 1, 18, 59, 44, 32, 52, 27, 12, 6, 17, 60, 45, 37, 49, 28, 19, 2, 15, 61, 47, 34, 48, 29, 16, 5, 10, 62, 50, 33, 46, 30, 21, 8, 9, 63, 51, 36, 43, 31, 20, 7, 4}},
       /* O3 */ {16, 4, {60, 53, 45, 35, 37, 42, 49, 58, 28, 23, 19, 13, 2, 0, 15, 25, 61, 56, 47, 38, 34, 40, 48, 55, 29, 24, 16, 11, 5, 3, 10, 22, 62, 57, 50, 41, 33, 39, 46, 54, 30, 26, 21, 14, 8, 1, 9, 18, 63, 59, 51, 44, 36, 32, 43, 52, 31, 27, 20, 12, 7, 6, 4, 17}},
       /* O4 */ {16, 4, {17, 4, 6, 7, 12, 20, 27, 31, 52, 43, 32, 36, 44, 51, 59, 63, 18, 9, 1, 8, 14, 21, 26, 30, 54, 46, 39, 33, 41, 50, 57, 62, 22, 10, 3, 5, 11, 16, 24, 29, 55, 48, 40, 34, 38, 47, 56, 61, 25, 15, 0, 2, 13, 19, 23, 28, 58, 49, 42, 37, 35, 45, 53, 60}},
       /* P3 */ {14, 5, {60, 53, 45, 35, 32, 46, 55, 28, 20, 14, 5, 0, 15, 25, 61, 56, 47, 38, 37, 43, 54, 29, 23, 12, 8, 3, 10, 22, 62, 57, 50, 41, 34, 42, 52, 30, 24, 19, 7, 1, 9, 18, 63, 59, 51, 44, 33, 40, 49, 31, 26, 16, 13, 6, 4, 17, -1, -1, -1, -1, 36, 39, 48, 58, 27, 21, 11, 2, -1, -1}},
       /* P4 */ {14, 5, {60, 53, 44, 33, 40, 49, 31, 26, 16, 13, 2, 0, 15, 25, 61, 56, 45, 36, 39, 48, 58, 27, 21, 11, 5, 3, 10, 22, 62, 57, 47, 35, 32, 46, 55, 28, 20, 14, 8, 1, 9, 18, 63, 59, 50, 38, 37, 43, 54, 29, 23, 12, 7, 6, 4, 17, -1, -1, 51, 41, 34, 42, 52, 30, 24, 19, -1, -1, -1, -1}},
       /* Q3 */ {16, 5, {-1, -1, 56, 45, 36, 39, 48, 58, 28, 23, 19, 13, 2, 0, 15, 25, -1, -1, 57, 47, 35, 32, 46, 55, 29, 24, 16, 11, 5, 3, 10, 22, -1, -1, 59, 50, 38, 37, 43, 54, 30, 26, 21, 14, 8, 1, 9, 18, -1, -1, 60, 51, 41, 34, 42, 52, 31, 27, 20, 12, 7, 6, 4, 17, 63, 62, 61, 53, 44, 33, 40, 49, -1, -1, -1, -1, -1, -1, -1, -1}},
       /* Q4 */ {16, 5, {60, 53, 45, 35, 37, 42, 49, 58, 27, 21, 11, 2, 4, 18, -1, -1, 61, 56, 47, 38, 34, 40, 48, 55, 28, 20, 14, 5, 0, 17, -1, -1, 62, 57, 50, 41, 33, 39, 46, 54, 29, 23, 12, 8, 3, 15, -1, -1, 63, 59, 51, 44, 36, 32, 43, 52, 30, 24, 19, 7, 1, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 31, 26, 16, 13, 6, 9, 22, 25}}},
      /* PS */
      {{0.714285714, 2.5}, {0.714285714, 5}}};
  }
}
class CathodeSegmentationCreatorRegisterCreateSegType2
{
 public:
  CathodeSegmentationCreatorRegisterCreateSegType2()
  {
    registerCathodeSegmentationCreator(2, createSegType2);
  }
} aCathodeSegmentationCreatorRegisterCreateSegType2;

} // namespace impl4
} // namespace mapping
} // namespace mch
} // namespace o2
