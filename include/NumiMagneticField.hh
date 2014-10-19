
// $Id: NumiMagneticField.hh,v 1.2.4.4 2014/10/19 00:15:58 lebrun Exp $
// --------------------------------------------------------------
// NumiMagneticField.hh modified by Yuki 2004/7/16
// modified by Yuki 8/2/04

#ifndef NumiMagneticField_H
#define NumiMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "NumiDataInput.hh"
class NumiMagneticField : public G4MagneticField
{
  public:
    NumiMagneticField();
    ~NumiMagneticField();

  virtual void GetFieldValue( const double Point[3],  double *Bfield ) const;  
  void DumpItToAsciiFile(); // Not implemented, but could be useful..
  
  private:
    NumiDataInput* NumiData;
    mutable bool dumpHasBeenDump;
  //G4double current;
};

class NumiMagneticFieldIC : public G4MagneticField
{
  public:
    NumiMagneticFieldIC();
    ~NumiMagneticFieldIC();

  virtual void GetFieldValue( const double Point[3], double *Bfield ) const;
  void DumpItToAsciiFile();

  private:
    NumiDataInput* NumiData;
    mutable bool dumpHasBeenDump;
  //G4double current;
};

class NumiMagneticFieldOC : public G4MagneticField
{
  public:
    NumiMagneticFieldOC();
    ~NumiMagneticFieldOC();

  virtual void GetFieldValue( const double Point[3], double *Bfield ) const;

  private:
    NumiDataInput* NumiData;
  //G4double current;
};

#endif

