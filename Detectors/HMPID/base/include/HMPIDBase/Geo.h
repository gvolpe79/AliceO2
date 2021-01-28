// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_HMPID_GEO_H
#define ALICEO2_HMPID_GEO_H

#include <TMath.h>

namespace o2
{
namespace hmpid
{

/*  ------------------ HMPID Detector Coordinate definition -------------------

  143 .-------------.             .-------------.          0  ------.------   23
      |             |             |.-----.-----.|            |      | -->  |
      |             |             ||     |     ||            |      | 0  9 |
      |             |             ||  4  |  5  ||            |      | Dil. |
    ^ |             |             ||     |     ||          | |      |      |  ^
    | |             |             ||_____|_____||          | |      |      |  |
    | |             |             |.-----.-----.|          | |      |      |  |
    Y |             |             ||     |     ||            |      |      |
    | |             |             ||  2  |  3  ||     Column |  n   |  n+1 | Column
    | |             |             ||     |     ||            |      |      |
      |             |             ||_____|_____||          | |      |      |  |
      |             |       47  ^ |.-----.-----.|          | |      |      |  |
      |             |           | ||     |     ||          V |      |      |  |
      |             |           Y ||  0  |  1  ||            | Dil. |      |
      |             |           | ||     |     ||            | 9  0 |      |
      |             |        0    ||_____|_____||            |  <-- |      |
    0 .-------------.              -------------          23  ------.------   0
      0  --X-->   159             0 -X->79

     Pad(Module,x,y)         Pad(Chamber,PhotoCat,x,y)  Pad(Equipment,Column,Dilogic,Channel)

            Equipment n              Equipment n+1

 143
  ^
  |    40 32 24 16 08 00           07 15 23 31 39 47
  |    41 33 25 17 09 01           06 14 22 30 38 46
       42 34 26 18 10 02           05 13 21 29 37 45
  y    43 35 27 19 11 03           04 12 20 28 36 44  (Column,Dil.)
       44 36 28 20 12 04           03 11 19 27 35 43
  |    45 37 29 21 13 05           02 10 18 26 34 42
  |    46 38 30 22 14 06           01 09 17 25 33 41
       47 39 31 23 15 07           00 08 16 24 32 40
  0


   0 --------- x ------> 79      80 ------ x -----------> 159

     For Equipment n :    x = 79 - (Dilo * 8 + Chan / 8)
                          y = 143 - (Column * 6 + Chan % 6)

     For Equipment n+1 :  x = 80 + (Dilo * 8 + Chan / 8)
                          y = Column * 6 + Chan % 6


    --------------------------------------------------------------------------- */
/// \class Geo
/// \brief HMPID  detector geometry (only statics)
class Geo
{
 public:
  // From AliTOFGeometry
//  static void translate(Float_t* xyz, Float_t translationVector[3]);
//  enum {
//    // DAQ characteristics
//    kNDDL = 4,    // Number of DDL (Detector Data Link) per sector
//    kNTRM = 12,   // Number of TRM ( Readout Module) per DDL
//    kNTdc = 15,   // Number of Tdc (Time to Digital Converter) per TRM
//    kNChain = 2,  // Number of chains per TRM
//    kNCrate = 72, // Number of Crates
//    kNCh = 8      // Number of channels per Tdc
//  };


  // ---- HMPID geometry -------
  static constexpr int MAXEQUIPMENTS = 14;
  static constexpr int N_SEGMENTS = 3;
  static constexpr int N_COLXSEGMENT = 8;
  static constexpr int N_COLUMNS = 24;
  static constexpr int N_DILOGICS = 10;
  static constexpr int N_CHANNELS = 48;
  static constexpr int N_DILOCHANNELS = 64;

  static constexpr int N_MODULES = 7;
  static constexpr int N_XROWS = 160;
  static constexpr int N_YCOLS = 144;

  static constexpr int MAXYCOLS = 143;
  static constexpr int MAXHALFXROWS = 79;
  static constexpr int HALFXROWS = 80;

  static constexpr int DILOPADSCOLS = 6;
  static constexpr int DILOPADSROWS = 8;

  static constexpr int EQUIPMENTSPERMODULE = 2;

  static constexpr int N_EQUIPMENTTOTALPADS = N_SEGMENTS * N_COLXSEGMENT * N_DILOGICS * N_CHANNELS;
  static constexpr int N_HMPIDTOTALPADS = MAXEQUIPMENTS * N_SEGMENTS * N_COLXSEGMENT * N_DILOGICS * N_CHANNELS;

  static constexpr int N_PHOTOCATODS = 6;
  static constexpr int N_PHOTOCATODSX = 80;
  static constexpr int N_PHOTOCATODSY = 48;
  static constexpr int MAXXPHOTO = 79;
  static constexpr int MAXYPHOTO = 47;


  static void Module2Equipment(int Mod, int Col, int Row, int *Equi, int *Colu, int *Dilo, int *Chan);
  static void Equipment2Module(int Equi, int Colu, int Dilo, int Chan, int *Mod, int *Col, int *Row);
  
  
  // from
  //static constexpr Bool_t FEAWITHMASKS[NSECTORS] =
  //  // TOF sectors with Nino masks: 0, 8, 9, 10, 16
  //  {kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE,
  //   kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE, kFALSE};
  //; // Selecting TOF sectors containing FEA cooling masks


//  static Float_t getCableLength(Int_t icrate, Int_t islot, Int_t ichain, Int_t itdc) { return CABLELENGTH[icrate][islot - 3][ichain][itdc / 3]; }

 private:
  static void Init();

  ClassDefNV(Geo, 1);
};
} // namespace hmpid
} // namespace o2

#endif
