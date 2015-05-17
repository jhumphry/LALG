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
      Register_Routine (T, Check_Rot'Access,
                        "Check real rot and rotg Givens rotation routines.");
      Register_Routine (T, Check_Rotm'Access,
                        "Check real rotm and rotmg Modified Givens rotation routines.");
      Register_Routine (T, Check_Swap'Access, "Check real swap routine.");
      Register_Routine (T, Check_Scal'Access, "Check real scal routine.");
      Register_Routine (T, Check_Copy'Access, "Check real copy routine.");
      Register_Routine (T, Check_Asum'Access, "Check real asum routine.");
   end Register_Tests;

   ------------
   -- Set_Up --
   ------------

   procedure Set_Up (T : in out Level1_Test) is
   begin
      null;
   end Set_Up;

   ---------------
   -- Check_Rot --
   ---------------

   procedure Check_Rot (T : in out Test_Cases.Test_Case'Class) is
      A : Concrete_Real_Matrix := Make(((6.0, 5.0, 0.0),
                                        (5.0, 1.0, 4.0),
                                       (0.0, 4.0, 3.0)));
      R1 : Real_Matrix_Vector := A.Row(1);
      R2 : Real_Matrix_Vector := A.Row(2);
      Expected_Result : Real_2D_Array := ((7.8102, 4.4813, 2.5607),
                                          (0.0000,-2.4327, 3.0729),
                                          (0.0000, 4.0000, 3.0000));
      a1, b1 : Float;
      c, s : Float;
   begin
      a1 := A(1,1); b1 := A(2,1);
      rotg(a1, b1, c, s);
      Assert(abs(a1 - 7.8102) <= 0.001, "Givens rotation r not set correctly");
      Assert(abs(c - 0.7682) <= 0.001, "Givens rotation c not set correctly");
      Assert(abs(s - 0.6402) <= 0.001, "Givens rotation s not set correctly");
      rot(R1, R2, c, s);
      Assert(Approx_Equal(A, Expected_Result, 0.001), "Givens rotation not applied correctly");

   end Check_Rot;

   ----------------
   -- Check_Rotm --
   ----------------

   procedure Check_Rotm (T : in out Test_Cases.Test_Case'Class) is
      A : Concrete_Real_Matrix := Make(((6.0, 5.0, 0.0),
                                        (5.0, 1.0, 4.0),
                                       (0.0, 4.0, 3.0)));
      R1 : Real_Matrix_Vector := A.Row(1);
      R2 : Real_Matrix_Vector := A.Row(2);
      Expected_Result : Real_2D_Array := ((7.8102, 4.4813, 2.5607),
                                          (0.0000,-2.4327, 3.0729),
                                          (0.0000, 4.0000, 3.0000));
      a1, b1, d1, d2 : Float;
      Params : Modified_Givens_Params;
   begin
      d1 := 0.707; d2 := 0.707;
      a1 := A(1,1); b1 := A(2,1);
      rotmg(d1, d2, a1, b1, Params);
      rotm(R1, R2, Params);
      Assert(abs(A(2,1))<0.001, "Modified Givens rotation not applied correctly");
   end Check_Rotm;

   ----------------
   -- Check_Swap --
   ----------------

   procedure Check_Swap (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      swap(X, Y);
      Assert(X = (4.0, 5.0, 6.0), "X not swapped correctly.");
      Assert(Y = (1.0, 2.0, 3.0), "Y not swapped correctly.");
      swap(X2, Y2);
      Assert(X = (4.0, 1.0, 3.0), "X not swapped correctly via view.");
      Assert(Y = (5.0, 2.0, 6.0), "Y not swapped correctly via view.");
   end Check_Swap;

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
   -- Check_Copy --
   ----------------

   procedure Check_Copy (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      copy(X, Y);
      Assert(X = (1.0, 2.0, 3.0), "X modified by a copy operation.");
      Assert(Y = (1.0, 2.0, 3.0), "Y not copied correctly.");
      copy(X2, Y2);
      Assert(X = (1.0, 2.0, 3.0), "X modified by a copy via view.");
      Assert(Y = (2.0, 2.0, 3.0), "Y not copied correctly via view.");
   end Check_Copy;

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
