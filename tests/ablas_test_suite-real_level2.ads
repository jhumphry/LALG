-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 2 BLAS routines

with AUnit; use Aunit;
with AUnit.Test_Cases; use AUnit.Test_Cases;

with aBLAS.Real_BLAS.Util;
private with Ada.Strings.Unbounded;

generic
   with package BLAS is new aBLAS.Real_BLAS(<>);
   Real_Type_Name : String := "<Unknown>";
   Epsilon : BLAS.Real := BLAS."*"(2.0, BLAS.Real'Model_Epsilon);
   Soft_Epsilon : BLAS.Real := BLAS."*"(1000.0, BLAS.Real'Model_Epsilon);
package aBLAS_Test_Suite.Real_Level2 is

   use BLAS;
   package BLAS_Util is new BLAS.Util;
   use BLAS_Util;

   type Level2_Test is new Test_Cases.Test_Case with null record;

   procedure Register_Tests (T: in out Level2_Test);

   function Name (T : Level2_Test) return Test_String is
     (Format("Test real Level 2 BLAS routines for " &
               Real_Type_Name));

   procedure Set_Up (T : in out Level2_Test);

   procedure Check_Gemv (T : in out Test_Cases.Test_Case'Class);

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
     ( 1 => (Check_Gemv'Access, +"Check real gemv routine.") );

end aBLAS_Test_Suite.Real_Level2;
