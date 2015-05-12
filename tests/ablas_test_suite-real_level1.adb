-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 1 BLAS routines

with AUnit.Assertions; use AUnit.Assertions;

with aBLAS_Test_Suite;
use aBLAS_Test_Suite.Float_BLAS;

package body aBLAS_Test_Suite.Real_Level1 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T: in out Level1_Test) is
      use AUnit.Test_Cases.Registration;
   begin
      Register_Routine (T, Check_Scal'Access, "Check real scal routine.");
      Register_Routine (T, Check_Asum'Access, "Check real asum routine.");
   end Register_Tests;

   ------------
   -- Set_Up --
   ------------

   procedure Set_Up (T : in out Level1_Test) is
   begin
      null;
   end Set_Up;

   ----------------
   -- Check_Scal --
   ----------------

   procedure Check_Scal (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
   begin
      Assert(X = (1.0, 2.0, 3.0), "Vector equality");
      Assert(Y = (2.0, 3.0), "Vector equality on views");
      scal(X, 1.0);
      Assert(X = (1.0, 2.0, 3.0), "Scaling by 1.0");
      Assert(Y = (2.0, 3.0), "Scaling by 1.0 with views");
      scal(X, -1.0);
      Assert(X = (-1.0, -2.0, -3.0), "Scaling by -1.0");
      scal(X, -2.0);
      Assert(X = (2.0, 4.0, 6.0), "Scaling by -2.0");
      scal(X, 100.0);
      Assert(X = (200.0, 400.0, 600.0), "Scaling by 100.0");
      Assert(Y = (400.0, 600.0), "Scaling by 100.0 on views");
      scal(Y, 0.25);
      Assert(X = (200.0, 100.0, 150.0), "Scaling by 0.25 via viewsk");
      Assert(Y = (100.0, 150.0), "Scaling by 0.25 on views");
   end Check_Scal;

   ----------------
   -- Check_Asum --
   ----------------

   procedure Check_Asum (T : in out Test_Cases.Test_Case'Class) is
      X : Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : Concrete_Real_Vector := Make((-4.0, 5.0, -6.0));
   begin
      Assert(asum(X) = 6.0, "asum(X) not working");
      Assert(asum(Y) = 15.0, "asum(Y) not working (not absolute values?)");
   end Check_Asum;

end aBLAS_Test_Suite.Real_Level1;
