-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 1 BLAS routines

with AUnit; use Aunit;
with AUnit.Test_Cases; use AUnit.Test_Cases;

package aBLAS_Test_Suite.Real_Level1 is

   type Level1_Test is new Test_Cases.Test_Case with null record;

   procedure Register_Tests (T: in out Level1_Test);

   function Name (T : Level1_Test) return Test_String is
      (Format("Test real Level 1 BLAS routines"));

   procedure Set_Up (T : in out Level1_Test);

   procedure Check_Rot (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Rotm (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Swap (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Scal (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Copy (T : in out Test_Cases.Test_Case'Class);
   procedure Check_Asum (T : in out Test_Cases.Test_Case'Class);

end aBLAS_Test_Suite.Real_Level1;
