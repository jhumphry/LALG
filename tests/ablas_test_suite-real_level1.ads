-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 1 BLAS routines

with AUnit; use Aunit;
with AUnit.Test_Cases; use AUnit.Test_Cases;

private with Ada.Strings.Unbounded;

generic
   with package BLAS is new aBLAS.Real_BLAS(<>);
   Real_Type_Name : String := "<Unknown>";
   Epsilon : BLAS.Real := BLAS."*"(2.0, BLAS.Real'Model_Epsilon);
   Soft_Epsilon : BLAS.Real := BLAS."*"(1000.0, BLAS.Real'Model_Epsilon);
package aBLAS_Test_Suite.Real_Level1 is

   use BLAS;

   type Level1_Test is new Test_Cases.Test_Case with null record;

   procedure Register_Tests (T: in out Level1_Test);

   function Name (T : Level1_Test) return Test_String is
     (Format("Test real Level 1 BLAS routines for " &
               Real_Type_Name));

   procedure Set_Up (T : in out Level1_Test);

   procedure Check_Rot (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Rotm (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Swap (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Scal (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Copy (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Axpy (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Dot_Sdsdot (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Nrm2 (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Asum (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Iamax (T : in out Test_Cases.Test_Case'Class);

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
     ( (Check_Rot'Access, +"Check real rot and rotg Givens rotation routines."),
       (Check_Rotm'Access, +"Check real rotm and rotmg Modified Givens rotation routines."),
       (Check_Swap'Access, +"Check real swap routine."),
       (Check_Scal'Access, +"Check real scal routine."),
       (Check_Copy'Access, +"Check real copy routine."),
       (Check_Axpy'Access, +"Check real axpy routine."),
       (Check_Dot_Sdsdot'Access, +"Check real dot and sdsdot routine."),
       (Check_Nrm2'Access, +"Check real nrm2 routine."),
       (Check_Asum'Access, +"Check real asum routine."),
       (Check_Iamax'Access, +"Check real iamax routine.") );

end aBLAS_Test_Suite.Real_Level1;
