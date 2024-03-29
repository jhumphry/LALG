-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines
-- Unit tests for Level 3 BLAS routines

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
package BLAS_Real_Level3 is

   pragma Unreferenced(Epsilon);
   pragma Unreferenced(Soft_Epsilon);

   use BLAS_Base, BLAS;

   type Level3_Test is new Test_Cases.Test_Case with null record;

   procedure Register_Tests (T: in out Level3_Test);

   function Name (T : Level3_Test) return Test_String is
     (Format("Test real Level 3 BLAS routines for " &
               Real_Type_Name));

   procedure Set_Up (T : in out Level3_Test);

   procedure Check_Gemm (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Symm (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Syrk (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Syr2k (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Trmm (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Trsm (T : in out Test_Cases.Test_Case'Class);

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
     ( (Check_Gemm'Access, +"Check real gemm routine."),
       (Check_Symm'Access, +"Check real symm routine."),
       (Check_Syrk'Access, +"Check real syrk routine."),
       (Check_Syr2k'Access, +"Check real syr2k routine."),
       (Check_Trmm'Access, +"Check real trmm routine."),
       (Check_Trsm'Access, +"Check real trsm routine.")
      );

end BLAS_Real_Level3;
