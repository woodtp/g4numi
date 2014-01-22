
// $Id: NumiMagneticField.hh,v 1.2.4.2 2014/01/22 22:31:06 kordosky Exp $
// --------------------------------------------------------------
// NumiMagneticField.hh modified by Yuki 2004/7/16
// modified by Yuki 8/2/04

#ifndef NumiMagneticField_H
#define NumiMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "NumiDataInput.hh"

#include "NumiKelvinFunctions.hh"

class NumiMagneticField : public G4MagneticField
{
  public:
    NumiMagneticField();
    ~NumiMagneticField();

  virtual void GetFieldValue( const double Point[3],
                               double *Bfield ) const;
  
  private:
    NumiDataInput* NumiData;
  //G4double current;
};

class NumiMagneticFieldIC : public G4MagneticField
{
  public:
    NumiMagneticFieldIC();
    ~NumiMagneticFieldIC();

  virtual void GetFieldValue( const double Point[3], double *Bfield ) const;

  private:
    NumiDataInput* NumiData;
    NumiKelvinFunctions* NumiKelvin;
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

