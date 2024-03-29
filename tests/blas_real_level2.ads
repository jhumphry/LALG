-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines
-- Unit tests for Level 2 BLAS routines

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

with AUnit; use AUnit;
with AUnit.Test_Cases; use AUnit.Test_Cases;

with LALG, LALG.Real_BLAS;

private with Ada.Strings.Unbounded;

generic
   with package BLAS_Base is new LALG(<>);
   with package BLAS is new BLAS_Base.Real_BLAS;
   Real_Type_Name : String := "<Unknown>";
   Epsilon : BLAS_Base.Real := BLAS_Base."*"(2.0, BLAS_Base.Real'Model_Epsilon);
   Soft_Epsilon : BLAS_Base.Real := BLAS_Base."*"(1000.0, BLAS_Base.Real'Model_Epsilon);
package BLAS_Real_Level2 is

   pragma Unreferenced(Epsilon);
   pragma Unreferenced(Soft_Epsilon);

   use BLAS_Base, BLAS;

   type Level2_Test is new Test_Cases.Test_Case with null record;

   procedure Register_Tests (T: in out Level2_Test);

   function Name (T : Level2_Test) return Test_String is
     (Format("Test real Level 2 BLAS routines for " &
               Real_Type_Name));

   procedure Set_Up (T : in out Level2_Test);

   procedure Check_Gemv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Symv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Spmv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Trmv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Tpmv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Trsv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Tpsv (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Ger (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Syr_Syr2 (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Spr_Spr2 (T : in out Test_Cases.Test_Case'Class);

private

   -- RM 3.10.2(32) means we cannot hide this away inside the body of the
   -- generic unit - the use of unbounded strings and the Test_Details type
   -- is just to make things less irritating.

   use Ada.Strings.Unbounded;

   function "+"(Source : in String) return Unbounded_String
                renames To_Unbounded_String;

   type Test_Details is
      record
         T : Test_Routine;
         D : Unbounded_String;
      end record;

   Test_Details_List: array (Positive range <>) of Test_Details :=
     ( (Check_Gemv'Access, +"Check real gemv routine."),
       (Check_Symv'Access, +"Check real symv routine."),
       (Check_Spmv'Access, +"Check real spmv routine."),
       (Check_Trmv'Access, +"Check real trmv routine."),
       (Check_Tpmv'Access, +"Check real tpmv routine."),
       (Check_Trsv'Access, +"Check real trsv routine."),
       (Check_Tpsv'Access, +"Check real tpsv routine."),
       (Check_Ger'Access, +"Check real ger routine."),
       (Check_Syr_Syr2'Access, +"Check real syr and syr2 routines."),
       (Check_Spr_Spr2'Access, +"Check real spr and spr2 routines.")
      );

end BLAS_Real_Level2;
