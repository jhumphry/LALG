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
      Register_Routine (T, Check_Axpy'Access, "Check real axpy routine.");
      Register_Routine (T, Check_Dot_Sdsdot'Access, "Check real dot and sdsdot routine.");
      Register_Routine (T, Check_Nrm2'Access, "Check real nrm2 routine.");
      Register_Routine (T, Check_Asum'Access, "Check real asum routine.");
      Register_Routine (T, Check_Iamax'Access, "Check real iamax routine.");
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
   -- Check_Axpy --
   ----------------

   procedure Check_Axpy (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      axpy(X, Y, 1.0);
      Assert(X = (1.0, 2.0, 3.0), "X modified by an Y<-aX+Y op.");
      Assert(Y = (5.0, 7.0, 9.0), "Y not correct following an Y<-aX+Y op.");
      axpy(X2, Y2, 2.0);
      Assert(X = (1.0, 2.0, 3.0), "X modified by an Y<-aX+Y op via view.");
      Assert(Y = (9.0, 7.0, 15.0), "Y not correct following an Y<-aX+Y op via view.");
   end Check_Axpy;

   ----------------------
   -- Check_Dot_Sdsdot --
   ----------------------

   procedure Check_Dot_Sdsdot (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      Assert(dot(X, Y) = 32.0, "dot incorrect.");
      Assert(sdsdot(X, Y, 4.0) = 36.0, "sdsdot incorrect.");
      Assert(dot(X2, Y2) = 26.0, "dot incorrect applied to views.");
      Assert(sdsdot(X2, Y2, -4.0) = 22.0, "sdsdot incorrect applied to views.");
   end Check_Dot_Sdsdot;

   ----------------
   -- Check_Nrm2 --
   ----------------

   procedure Check_Nrm2 (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
   begin
      Assert(abs(nrm2(X)-3.7416573867739413)<0.000001, "Nrm2 incorrect.");
      Assert(abs(nrm2(X2)-3.605551275463989)<0.000001, "Nrm2 incorrect applied to view.");
   end Check_Nrm2;

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

   -----------------
   -- Check_Iamax --
   -----------------

   procedure Check_Iamax (T : in out Test_Cases.Test_Case'Class) is
      X : Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((4.0, -5.0, -6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 1, Length => 2);
   begin
      Assert(iamax(X) = 3, "iamax(X) not working");
      Assert(iamax(Y) = 3, "iamax(Y) not working");
      Assert(iamax(Y2) = 2, "iamax(Y2) not working on view");
   end Check_Iamax;

end aBLAS_Test_Suite.Real_Level1;
