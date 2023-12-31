/*
 * Simple example c program to write a
 * binary datafile for tecplot.  This example
 * does the following:
 *
 *   1.  Open a datafile called "t.plt"
 *   2.  Assign values for X,Y, and P
 *   3.  Write out a zone dimensioned 4x5
 *   4.  Close the datafile.
 */

#include "TECIO.h"

#ifndef NULL
#define NULL 0
#endif

enum FileType { FULL = 0, GRID = 1, SOLUTION = 2 };

main ()
{
  float X[5][4], Y[5][4], P[5][4];
  double SolTime;
  INTEGER4 Debug,I,J,III,DIsDouble,VIsDouble,IMax,JMax,KMax,ZoneType,StrandID,ParentZn,IsBlock;
  INTEGER4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn, FileType;

  Debug     = 1;
  VIsDouble = 0;
  DIsDouble = 0;
  IMax      = 4;
  JMax      = 5;
  KMax      = 1;
  ZoneType  = 0;      /* Ordered */
  SolTime   = 360.0;
  StrandID  = 0;     /* StaticZone */
  ParentZn  = 0;      /* No Parent */
  IsBlock   = 1;      /* Block */
  ICellMax  = 0;
  JCellMax  = 0;
  KCellMax  = 0;
  NFConns   = 0;
  FNMode    = 0;
  ShrConn   = 0;
  FileType  = FULL;

/*
 * Open the file and write the tecplot datafile 
 * header information 
 */
  I = TECINI111("SIMPLE DATASET",
                "X Y P",
                "t.plt",
                ".",
                &FileType,
                &Debug,
                &VIsDouble);

  for (J = 0; J < 5; J++)
  for (I = 0; I < 4; I++)
    {
      X[J][I] = (float)(I+1);
      Y[J][I] = (float)(J+1);
      P[J][I] = (float)((I+1)*(J+1));
    }
/*
 * Write the zone header information.
 */
  I = TECZNE111("Simple Zone",
                &ZoneType,
                &IMax,
                &JMax,
                &KMax,
                &ICellMax,
                &JCellMax,
                &KCellMax,
                &SolTime,
                &StrandID,
                &ParentZn,
                &IsBlock,
                &NFConns,
                &FNMode,
                0,              /* TotalNumFaceNodes */
                0,              /* NumConnectedBoundaryFaces */
                0,              /* TotalNumBoundaryConnections */
                NULL,           /* PassiveVarList */
                NULL,           /* ValueLocation = Nodal */
                NULL,           /* SharVarFromZone */
                &ShrConn);
/*
 * Write out the field data.
 */
  III = IMax*JMax;
  I   = TECDAT111(&III,&X[0][0],&DIsDouble);
  I   = TECDAT111(&III,&Y[0][0],&DIsDouble);
  I   = TECDAT111(&III,&P[0][0],&DIsDouble);

  I = TECEND111();
}
