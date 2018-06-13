// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "../../base/include/HMPIDBase/Param.h"
//#include "HMPIDBase/Param.h"
//#include "TGeoMatrix.h" 
#include <TLatex.h>         //TestTrans()  
#include <TView.h>          //TestTrans()
#include <TPolyMarker3D.h>  //TestTrans()
#include <TRotation.h>
#include <TParticle.h>      //Stack()    
#include <TGeoPhysicalNode.h> //ctor
#include <TGeoBBox.h>
#include <TF1.h>                 //ctor
#include <iostream>

using namespace o2::hmpid;

ClassImp(o2::hmpid::Param);

// Mathieson constant definition
const Double_t Param::fgkD     = 0.222500;  // ANODE-CATHODE distance 0.445/2
//                                                                                          K3 = 0.66 along the wires (anode-cathode/wire pitch=0.5625)
const Double_t Param::fgkSqrtK3x = TMath::Sqrt(0.66);
const Double_t Param::fgkK2x     = TMath::PiOver2()*(1 - 0.5*fgkSqrtK3x);
const Double_t Param::fgkK1x     = 0.25*fgkK2x*fgkSqrtK3x/TMath::ATan(fgkSqrtK3x);
const Double_t Param::fgkK4x     = fgkK1x/(fgkK2x*fgkSqrtK3x);
//                                                                                          K3 = 0.87 along the wires (anode-cathode/wire pitch=0.5625)
const Double_t Param::fgkSqrtK3y = TMath::Sqrt(0.87);
const Double_t Param::fgkK2y     = TMath::PiOver2()*(1 - 0.5*fgkSqrtK3y);
const Double_t Param::fgkK1y     = 0.25*fgkK2y*fgkSqrtK3y/TMath::ATan(fgkSqrtK3y);
const Double_t Param::fgkK4y     = fgkK1y/(fgkK2y*fgkSqrtK3y);
//
  

Float_t Param::fgkMinPcX[]={0.,0.,0.,0.,0.,0.};
Float_t Param::fgkMaxPcX[]={0.,0.,0.,0.,0.,0.};
Float_t Param::fgkMinPcY[]={0.,0.,0.,0.,0.,0.};
Float_t Param::fgkMaxPcY[]={0.,0.,0.,0.,0.,0.};

Bool_t Param::fgMapPad[160][144][7];

Float_t Param::fgCellX=0.;
Float_t Param::fgCellY=0.;

Float_t Param::fgPcX=0;
Float_t Param::fgPcY=0;

Float_t Param::fgAllX=0;
Float_t Param::fgAllY=0;

Bool_t Param::fgInstanceType=kTRUE;  

Param* Param::fgInstance=0x0;        //singleton pointer               

Int_t Param::fgNSigmas  = 4;
Int_t Param::fgThreshold= 4;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Param::Param(Bool_t noGeo):
  fX(0), fY(0), fRefIdx(1.28947),fPhotEMean(6.675),fTemp(25)                          //just set a refractive index for C6F14 at ephot=6.675 eV @ T=25 C
{
// Here all the intitializition is taken place when Param::Instance() is invoked for the first time.
// In particular, matrices to be used for LORS<->MARS trasnformations are initialized from TGeo structure.    
// Note that TGeoManager should be already initialized from geometry.root file  

 /* AliCDBManager *pCDB = AliCDBManager::Instance();
  if(!pCDB) {
     AliWarning("No Nmean C6F14 from OCDB. Default is taken from ctor.");
  } else {
    AliCDBEntry *pNmeanEnt =pCDB->Get("HMPID/Calib/Nmean"); //contains TObjArray of 42 TF1 + 1 EPhotMean
    if(!pNmeanEnt) {
      AliWarning("No Nmean C6F14 from OCDB. Default is taken from ctor.");
    } else {
      TObjArray *pNmean = (TObjArray*)pNmeanEnt->GetObject();
      if(pNmean->GetEntries()==43) {                                               //for backward compatibility
        Double_t tmin,tmax;
        ((TF1*)pNmean->At(42))->GetRange(tmin,tmax);
        fPhotEMean = ((TF1*)pNmean->At(42))->Eval(tmin);                          //photon eMean from OCDB
        AliInfo(Form("EPhotMean = %f eV successfully loaded from OCDB",fPhotEMean));
      } else {
        AliWarning("For backward compatibility EPhotMean is taken from ctor.");
      }
    }
  }
*/
  fRefIdx = MeanIdxRad(); //initialization of the running ref. index of freon
  
  Float_t dead=2.6;// cm of the dead zones between PCs-> See 2CRC2099P1


  if(noGeo==kTRUE) fgInstanceType=kFALSE;                                                   //instance from ideal geometry, no actual geom is present
    
  if(noGeo==kFALSE && !gGeoManager)  
  {
    TGeoManager::Import("geometry.root");
    if(!gGeoManager) Printf("!!!!!!No geometry loaded!!!!!!!");
  }
  
  fgCellX=0.8;fgCellY=0.84;
  
  if(!noGeo==kTRUE){
    TGeoVolume *pCellVol = gGeoManager->GetVolume("Hcel");
    if(pCellVol) {
      TGeoBBox *bcell = (TGeoBBox *)pCellVol->GetShape();
      fgCellX=2.*bcell->GetDX(); fgCellY = 2.*bcell->GetDY();  // overwrite the values with the read ones
    }
  }    
  fgPcX=80.*fgCellX; fgPcY = 48.*fgCellY;
  fgAllX=2.*fgPcX+dead;
  fgAllY=3.*fgPcY+2.*dead;

  fgkMinPcX[1]=fgPcX+dead; fgkMinPcX[3]=fgkMinPcX[1];  fgkMinPcX[5]=fgkMinPcX[3];
  fgkMaxPcX[0]=fgPcX; fgkMaxPcX[2]=fgkMaxPcX[0];  fgkMaxPcX[4]=fgkMaxPcX[2];
  fgkMaxPcX[1]=fgAllX; fgkMaxPcX[3]=fgkMaxPcX[1];  fgkMaxPcX[5]=fgkMaxPcX[3];

  fgkMinPcY[2]=fgPcY+dead; fgkMinPcY[3]=fgkMinPcY[2];  
  fgkMinPcY[4]=2.*fgPcY+2.*dead; fgkMinPcY[5]=fgkMinPcY[4];
  fgkMaxPcY[0]=fgPcY; fgkMaxPcY[1]=fgkMaxPcY[0];  
  fgkMaxPcY[2]=2.*fgPcY+dead; fgkMaxPcY[3]=fgkMaxPcY[2]; 
  fgkMaxPcY[4]=fgAllY; fgkMaxPcY[5]=fgkMaxPcY[4];   
    
  fX=0.5*SizeAllX();
  fY=0.5*SizeAllY();
  
      
  for(Int_t ich=kMinCh;ich<=kMaxCh;ich++) {
    for(Int_t padx=0;padx<160;padx++) {
       for(Int_t pady=0;pady<144;pady++) {
         fgMapPad[padx][pady][ich] = kTRUE;             //init all the pads are active at the beginning....
       }
     }
   }
     

  for(Int_t i=kMinCh;i<=kMaxCh;i++)
    if(gGeoManager && gGeoManager->IsClosed()) {
      TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(Form("/HMPID/Chamber%i",i));
      if (!pne) {
        //AliErrorClass(Form("The symbolic volume %s does not correspond to any physical entry!",Form("HMPID_%i",i)));
        fM[i]=new TGeoHMatrix;
        IdealPosition(i,fM[i]);
      } else {
        TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
        if(pnode) fM[i]=new TGeoHMatrix(*(pnode->GetMatrix()));
        else {
          fM[i]=new TGeoHMatrix;
          IdealPosition(i,fM[i]);
        }
      }
    } else{
      fM[i]=new TGeoHMatrix;
      IdealPosition(i,fM[i]);
    } 
  fgInstance=this; 
}//ctor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::Print(Option_t* opt) const
{
// print some usefull (hopefully) info on some internal guts of HMPID parametrisation 
  
  for(Int_t i=0;i<7;i++) fM[i]->Print(opt);
}//Print()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::IdealPosition(Int_t iCh, TGeoHMatrix *pMatrix)
{
// Construct ideal position matrix for a given chamber
// Arguments: iCh- chamber ID; pMatrix- pointer to precreated unity matrix where to store the results
//   Returns: none
  const Double_t kAngHor=19.5;        //  horizontal angle between chambers  19.5 grad
  const Double_t kAngVer=20;          //  vertical angle between chambers    20   grad     
  const Double_t kAngCom=30;          //  common HMPID rotation with respect to x axis  30   grad     
  const Double_t kTrans[3]={490,0,0}; //  center of the chamber is on window-gap surface
  pMatrix->RotateY(90);               //  rotate around y since initial position is in XY plane -> now in YZ plane
  pMatrix->SetTranslation(kTrans);    //  now plane in YZ is shifted along x 
  switch(iCh){
    case 0:                pMatrix->RotateY(kAngHor);  pMatrix->RotateZ(-kAngVer);  break; //right and down 
    case 1:                                            pMatrix->RotateZ(-kAngVer);  break; //down              
    case 2:                pMatrix->RotateY(kAngHor);                               break; //right 
    case 3:                                                                         break; //no rotation
    case 4:                pMatrix->RotateY(-kAngHor);                              break; //left   
    case 5:                                            pMatrix->RotateZ(kAngVer);   break; //up
    case 6:                pMatrix->RotateY(-kAngHor); pMatrix->RotateZ(kAngVer);   break; //left and up 
  }
  pMatrix->RotateZ(kAngCom);     //apply common rotation  in XY plane    
   
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*Int_t Param::Stack(Int_t evt,Int_t tid)
{
// Prints some useful info from stack
// Arguments: evt - event number. if not -1 print info only for that event
//            tid - track id. if not -1 then print it and all it's mothers if any   
//   Returns: mother tid of the given tid if any
  AliRunLoader *pAL=AliRunLoader::Open(); 
  if(pAL->LoadHeader()) return -1;
  if(pAL->LoadKinematics()) return -1;
  
  Int_t mtid=-1;
  Int_t iNevt=pAL->GetNumberOfEvents();
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(evt!=-1 && evt!=iEvt) continue; //in case one needs to print the requested event, ignore all others
    pAL->GetEvent(iEvt);    
    AliStack *pStack=pAL->Stack();  
    if(tid==-1){                        //print all tids for this event
      for(Int_t i=0;i<pStack->GetNtrack();i++) pStack->Particle(i)->Print();
          Printf("totally %i tracks including %i primaries for event %i out of %i event(s)",
          pStack->GetNtrack(),pStack->GetNprimary(),iEvt,iNevt);
    }else{                              //print only this tid and it;s mothers
      if(tid<0 || tid>pStack->GetNtrack()) {Printf("Wrong tid, valid tid range for event %i is 0-%i",iEvt,pStack->GetNtrack());break;}
      TParticle *pTrack=pStack->Particle(tid); mtid=pTrack->GetFirstMother();
      TString str=pTrack->GetName();
      while((tid=pTrack->GetFirstMother()) >= 0){
        pTrack=pStack->Particle(tid);
        str+=" from ";str+=pTrack->GetName();
      } 
    }//if(tid==-1)      
  }//events loop
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return mtid;
}*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*Int_t Param::StackCount(Int_t pid,Int_t evt)
{
// Counts total number of particles of given sort (including secondary) for a given event
  AliRunLoader *pAL=AliRunLoader::Open(); 
  pAL->GetEvent(evt);    
  if(pAL->LoadHeader()) return 0;
  if(pAL->LoadKinematics()) return 0;
  AliStack *pStack=pAL->Stack();
  
  Int_t iCnt=0;
  for(Int_t i=0;i<pStack->GetNtrack();i++) if(pStack->Particle(i)->GetPdgCode()==pid) iCnt++;
  
  pAL->UnloadHeader();  pAL->UnloadKinematics();
  return iCnt;
}*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::Sigma2(Double_t trkTheta,Double_t trkPhi,Double_t ckovTh, Double_t ckovPh)
{
// Analithical calculation of total error (as a sum of localization, geometrical and chromatic errors) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  TVector3 v(-999,-999,-999);
  Double_t trkBeta = 1./(TMath::Cos(ckovTh)*GetRefIdx());
  
  if(trkBeta > 1) trkBeta = 1;                 //protection against bad measured thetaCer  
  if(trkBeta < 0) trkBeta = 0.0001;            //

  v.SetX(SigLoc (trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));
  v.SetY(SigGeom(trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));
  v.SetZ(SigCrom(trkTheta,trkPhi,ckovTh,ckovPh,trkBeta));

  return v.Mag2();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::SigLoc(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of localization error (due to finite segmentation of PC) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosf     = TMath::Cos(trkPhi);
  Double_t sinfd    = TMath::Sin(phiDelta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);        // formula (after 8 in the text)
  if (k<0) return 1e10;
  Double_t mu =sint*sinf+tantheta*(cost*cosfd*sinf+sinfd*cosf);                             // formula (10)
  Double_t e  =sint*cosf+tantheta*(cost*cosfd*cosf-sinfd*sinf);                             // formula (9)

  Double_t kk = betaM*TMath::Sqrt(k)/(GapThick()*alpha);                            // formula (6) and (7)
  Double_t dtdxc = kk*(k*(cosfd*cosf-cost*sinfd*sinf)-(alpha*mu/(betaM*betaM))*sint*sinfd); // formula (6)           
  Double_t dtdyc = kk*(k*(cosfd*sinf+cost*sinfd*cosf)+(alpha* e/(betaM*betaM))*sint*sinfd); // formula (7)            pag.4

  Double_t errX = 0.2,errY=0.25;                                                            //end of page 7
  return  TMath::Sqrt(errX*errX*dtdxc*dtdxc + errY*errY*dtdyc*dtdyc);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::SigCrom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of chromatic error (due to lack of knowledge of Cerenkov photon energy) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Fromulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    
  
  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                 // formula (11)
  Double_t dtdn = cost*GetRefIdx()*betaM*betaM/(alpha*tantheta);                    // formula (12)
            
//  Double_t f = 0.00928*(7.75-5.635)/TMath::Sqrt(12.);
  Double_t f = 0.0172*(7.75-5.635)/TMath::Sqrt(24.);

  return f*dtdn;
}//SigCrom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::SigGeom(Double_t trkTheta,Double_t trkPhi,Double_t thetaC, Double_t phiC,Double_t betaM)
{
// Analitical calculation of geometric error (due to lack of knowledge of creation point in radiator) on Cerenkov angle for a given Cerenkov photon 
// created by a given MIP. Formulae according to CERN-EP-2000-058 
// Arguments: Cerenkov and azimuthal angles for Cerenkov photon, [radians]
//            dip and azimuthal angles for MIP taken at the entrance to radiator, [radians]        
//            MIP beta
//   Returns: absolute error on Cerenkov angle, [radians]    

  Double_t phiDelta = phiC - trkPhi;

  Double_t sint     = TMath::Sin(trkTheta);
  Double_t cost     = TMath::Cos(trkTheta);
  Double_t sinf     = TMath::Sin(trkPhi);
  Double_t cosfd    = TMath::Cos(phiDelta);
  Double_t costheta = TMath::Cos(thetaC);
  Double_t tantheta = TMath::Tan(thetaC);
  
  Double_t alpha =cost-tantheta*cosfd*sint;                                                  // formula (11)
  
  Double_t k = 1.-GetRefIdx()*GetRefIdx()+alpha*alpha/(betaM*betaM);         // formula (after 8 in the text)
  if (k<0) return 1e10;

  Double_t eTr = 0.5*RadThick()*betaM*TMath::Sqrt(k)/(GapThick()*alpha);     // formula (14)
  Double_t lambda = (1.-sint*sinf)*(1.+sint*sinf);                                                  // formula (15)

  Double_t c1 = 1./(1.+ eTr*k/(alpha*alpha*costheta*costheta));                              // formula (13.a)
  Double_t c2 = betaM*TMath::Power(k,1.5)*tantheta*lambda/(GapThick()*alpha*alpha);  // formula (13.b)
  Double_t c3 = (1.+eTr*k*betaM*betaM)/((1+eTr)*alpha*alpha);                                // formula (13.c)
  Double_t c4 = TMath::Sqrt(k)*tantheta*(1-lambda)/(GapThick()*betaM);               // formula (13.d)
  Double_t dtdT = c1 * (c2+c3*c4);
  Double_t trErr = RadThick()/(TMath::Sqrt(12.)*cost);

  return trErr*dtdT;
}//SigGeom()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::SigmaCorrFact  (Int_t iPart, Double_t occupancy) 
{
  Double_t corr = 1.0;
                                                                                                             
  switch(iPart) {
    case 0: corr = 0.115*occupancy + 1.166; break; 
    case 1: corr = 0.115*occupancy + 1.166; break;
    case 2: corr = 0.115*occupancy + 1.166; break;
    case 3: corr = 0.065*occupancy + 1.137; break;
    case 4: corr = 0.048*occupancy + 1.202; break;
  }
                                                                                                                           
 return corr; 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Param* Param::Instance()
{
// Return pointer to the AliHMPIDParam singleton. 
// Arguments: none
//   Returns: pointer to the instance of AliHMPIDParam or 0 if no geometry       
  if(!fgInstance) new Param(kFALSE);                                //default setting for reconstruction, if no geometry.root -> AliFatal
  return fgInstance;  
}//Instance()    
 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Param* Param::InstanceNoGeo()
{
// Return pointer to the AliHMPIDParam singleton without the geometry.root. 
// Arguments: none
//   Returns: pointer to the instance of AliHMPIDParam or 0 if no geometry       
  if(!fgInstance) new Param(kTRUE);                               //to avoid AliFatal, for MOOD and displays, use ideal geometry parameters
  return fgInstance;  
}//Instance()    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t Param::IsInDead(Float_t x,Float_t y)
{
// Check is the current point is outside of sensitive area or in dead zones
// Arguments: x,y -position
//   Returns: 1 if not in sensitive zone           
  for(Int_t iPc=0;iPc<6;iPc++)
    if(x>=fgkMinPcX[iPc] && x<=fgkMaxPcX[iPc] && y>=fgkMinPcY[iPc] && y<=fgkMaxPcY [iPc]) return kFALSE; //in current pc
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t Param::IsDeadPad(Int_t padx,Int_t pady,Int_t ch)
{
// Check is the current pad is active or not
// Arguments: padx,pady pad integer coord
//   Returns: kTRUE if dead, kFALSE if active

    if(fgMapPad[padx-1][pady-1][ch]) return kFALSE; //current pad active
  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::Lors2Pad(Float_t x,Float_t y,Int_t &pc,Int_t &px,Int_t &py)
{
// Check the pad of given position
// Arguments: x,y- position [cm] in LORS; pc,px,py- pad where to store the result
//   Returns: none
  pc=px=py=-1;
  if     (x>fgkMinPcX[0] && x<fgkMaxPcX[0]) {pc=0; px=Int_t( x               / SizePadX());}//PC 0 or 2 or 4
  else if(x>fgkMinPcX[1] && x<fgkMaxPcX[1]) {pc=1; px=Int_t((x-fgkMinPcX[1]) / SizePadX());}//PC 1 or 3 or 5
  else return;
  if     (y>fgkMinPcY[0] && y<fgkMaxPcY[0]) {      py=Int_t( y               / SizePadY());}//PC 0 or 1
  else if(y>fgkMinPcY[2] && y<fgkMaxPcY[2]) {pc+=2;py=Int_t((y-fgkMinPcY[2]) / SizePadY());}//PC 2 or 3
  else if(y>fgkMinPcY[4] && y<fgkMaxPcY[4]) {pc+=4;py=Int_t((y-fgkMinPcY[4]) / SizePadY());}//PC 4 or 5
  else return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t Param::InHVSector(Float_t y)
{
//Calculate the HV sector corresponding to the cluster position
//Arguments: y
//Returns the HV sector in the single module
 
   Int_t hvsec = -1;
   Int_t pc,px,py;
   Lors2Pad(1.,y,pc,px,py);
   if(py==-1) return hvsec;
   
   hvsec = (py+(pc/2)*(kMaxPy+1))/((kMaxPy+1)/2);
   
   return hvsec;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t Param::FindTemp(Double_t tLow,Double_t tHigh,Double_t y)
{
//  Model for gradient in temperature
  Double_t yRad = HinRad(y);     //height in a given radiator
  if(tHigh<tLow) tHigh = tLow;   //if Tout < Tin consider just Tin as reference...
  if(yRad<0        ) yRad = 0;         //protection against fake y values
  if(yRad>SizePcY()) yRad = SizePcY(); //protection against fake y values
  
  Double_t gradT = (tHigh-tLow)/SizePcY();  // linear gradient
  return gradT*yRad+tLow;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::SetChStatus(Int_t ch,Bool_t status)
{
//Set a chamber on or off depending on the status
//Arguments: ch=chamber,status=kTRUE = active, kFALSE=off
//Returns: none
  for(Int_t padx=0;padx<kMaxPcx+1;padx++) {
     for(Int_t pady=0;pady<kMaxPcy+1;pady++) {
       fgMapPad[padx][pady][ch] = status;
     }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::SetSectStatus(Int_t ch,Int_t sect,Bool_t status)
{
//Set a given sector sect for a chamber ch on or off depending on the status
//Sector=0,5 (6 sectors)
//Arguments: ch=chamber,sect=sector,status: kTRUE = active, kFALSE=off
//Returns: none
  
  Int_t npadsect = (kMaxPcy+1)/6;
  Int_t padSectMin = npadsect*sect;
  Int_t padSectMax = padSectMin+npadsect;
  
  for(Int_t padx=0;padx<kMaxPcx+1;padx++) {
     for(Int_t pady=padSectMin;pady<padSectMax;pady++) {
       fgMapPad[padx][pady][ch] = status;
     }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::SetPcStatus(Int_t ch,Int_t pc,Bool_t status)
{
//Set a given PC pc for a chamber ch on or off depending on the status
//Arguments: ch=chamber,pc=PC,status: kTRUE = active, kFALSE=off
//Returns: none
  
  Int_t deltaX = pc%2;
  Int_t deltaY = pc/2;
  Int_t padPcXMin = deltaX*kPadPcX;
  Int_t padPcXMax = padPcXMin+kPadPcX;
  Int_t padPcYMin = deltaY*kPadPcY;
  Int_t padPcYMax = padPcYMin+kPadPcY;
  
  for(Int_t padx=padPcXMin;padx<padPcXMax;padx++) {
     for(Int_t pady=padPcYMin;pady<padPcYMax;pady++) {
       fgMapPad[padx][pady][ch] = status;
     }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::PrintChStatus(Int_t ch)
{
//Print the map status of a chamber on or off depending on the status
//Arguments: ch=chamber
//Returns: none
  Printf(" ");
  Printf(" --------- C H A M B E R  %d   ---------------",ch);
  for(Int_t pady=kMaxPcy;pady>=0;pady--) {
     for(Int_t padx=0;padx<kMaxPcx+1;padx++) {
       if(padx==80) printf(" ");
       printf("%d",fgMapPad[padx][pady][ch]);
     }
     printf(" %d \n",pady+1);
     if(pady%48==0) printf("\n");
   }
   printf("\n");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Param::SetGeomAccept()
{
//Set the real acceptance of the modules, due to ineficciency or hardware problems (up tp 1/6/2010)
//Arguments: none
//Returns: none
  SetSectStatus(0,3,kFALSE);
  SetSectStatus(4,0,kFALSE);
  SetSectStatus(5,1,kFALSE);
  SetSectStatus(6,2,kFALSE);
  SetSectStatus(6,3,kFALSE);
}
