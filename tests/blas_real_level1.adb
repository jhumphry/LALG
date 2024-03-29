-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines
-- Unit tests for Level 1 BLAS routines

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

with AUnit.Assertions; use AUnit.Assertions;

package body BLAS_Real_Level1 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T : in out Level1_Test) is
      use AUnit.Test_Cases.Registration;
   begin

      for I of Test_Details_List loop
         Register_Routine(T, I.T, To_String(I.D));
      end loop;

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
      pragma Unreferenced (T);
      A : constant Concrete_Real_Matrix := Make((
                                       (6.0, 5.0, 0.0),
                                       (5.0, 1.0, 4.0),
                                       (0.0, 4.0, 3.0)
                                      ));
      a1, b1 : Real;
      c, s : Real;
   begin
      a1 := A(1,1); b1 := A(2,1);
      rotg(a1, b1, c, s);
      Assert(abs(a1 - 7.810249675906654) <= Soft_Epsilon, "Givens rotation r not set correctly");
      Assert(abs(c - 0.76822127959737585) <= Soft_Epsilon, "Givens rotation c not set correctly");
      Assert(abs(s - 0.64018439966447993) <= Soft_Epsilon, "Givens rotation s not set correctly");

      -- TODO : test rot routine once high-precision result is available

      -- rot(R1, R2, c, s);
      -- Assert(Approx_Equal(A, Expected_Result, Soft_Epsilon), "Givens rotation not applied correctly");
   end Check_Rot;

   ----------------
   -- Check_Rotm --
   ----------------

   procedure Check_Rotm (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      A : Concrete_Real_Matrix := Make(((6.0, 5.0, 0.0),
                                        (5.0, 1.0, 4.0),
                                       (0.0, 4.0, 3.0)));
      R1 : Real_Matrix_Vector := A.Row(1);
      R2 : Real_Matrix_Vector := A.Row(2);
      a1, b1, d1, d2 : Real;
      Params : Modified_Givens_Params;
   begin
      d1 := 0.707; d2 := 0.707;
      a1 := A(1,1); b1 := A(2,1);
      rotmg(d1, d2, a1, b1, Params);
      rotm(R1, R2, Params);
      Assert(abs(A(2,1))<Epsilon, "Modified Givens rotation not applied correctly");
   end Check_Rotm;

   ----------------
   -- Check_Swap --
   ----------------

   procedure Check_Swap (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      swap(X, Y);
      Assert(X = Real_1D_Array'(4.0, 5.0, 6.0), "X not swapped correctly.");
      Assert(Y = Real_1D_Array'(1.0, 2.0, 3.0), "Y not swapped correctly.");
      swap(X2, Y2);
      Assert(X = Real_1D_Array'(4.0, 1.0, 3.0), "X not swapped correctly via view.");
      Assert(Y = Real_1D_Array'(5.0, 2.0, 6.0), "Y not swapped correctly via view.");
   end Check_Swap;

   ----------------
   -- Check_Scal --
   ----------------

   procedure Check_Scal (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
   begin
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "Vector equality");
      Assert(Y = Real_1D_Array'(2.0, 3.0), "Vector equality on views");
      scal(X, 1.0);
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "Scaling by 1.0");
      Assert(Y = Real_1D_Array'(2.0, 3.0), "Scaling by 1.0 with views");
      scal(X, -1.0);
      Assert(X = Real_1D_Array'(-1.0, -2.0, -3.0), "Scaling by -1.0");
      scal(X, -2.0);
      Assert(X = Real_1D_Array'(2.0, 4.0, 6.0), "Scaling by -2.0");
      scal(X, 100.0);
      Assert(X = Real_1D_Array'(200.0, 400.0, 600.0), "Scaling by 100.0");
      Assert(Y = Real_1D_Array'(400.0, 600.0), "Scaling by 100.0 on views");
      scal(Y, 0.25);
      Assert(X = Real_1D_Array'(200.0, 100.0, 150.0), "Scaling by 0.25 via viewsk");
      Assert(Y = Real_1D_Array'(100.0, 150.0), "Scaling by 0.25 on views");
   end Check_Scal;

   ----------------
   -- Check_Copy --
   ----------------

   procedure Check_Copy (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : constant Real_Vector_View := Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View := Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      copy(X, Y);
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X modified by a copy operation.");
      Assert(Y = Real_1D_Array'(1.0, 2.0, 3.0), "Y not copied correctly.");
      copy(X2, Y2);
      pragma Unreferenced(Y2);
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X modified by a copy via view.");
      Assert(Y = Real_1D_Array'(2.0, 2.0, 3.0), "Y not copied correctly via view.");
   end Check_Copy;

   ----------------
   -- Check_Axpy --
   ----------------

   procedure Check_Axpy (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : constant Real_Vector_View :=
        Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : Real_Vector_View :=
        Make(Y'Access, Start => 1, Stride => 2, Length => 2);
   begin
      axpy(X, Y, 1.0);
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X modified by an Y<-aX+Y op.");
      Assert(Y = Real_1D_Array'(5.0, 7.0, 9.0), "Y not correct following an Y<-aX+Y op.");
      axpy(X2, Y2, 2.0);
      Assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X modified by an Y<-aX+Y op via view.");
      Assert(Y = Real_1D_Array'(9.0, 7.0, 15.0), "Y not correct following an Y<-aX+Y op via view.");
   end Check_Axpy;

   ----------------------
   -- Check_Dot_Sdsdot --
   ----------------------

   procedure Check_Dot_Sdsdot (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      X2 : constant Real_Vector_View :=
        Make(X'Access, Start => 2, Stride => 1, Length => 2);
      Y : aliased Concrete_Real_Vector := Make((4.0, 5.0, 6.0));
      Y2 : constant Real_Vector_View :=
        Make(Y'Access, Start => 1, Stride => 2, Length => 2);
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
      pragma Unreferenced (T);
      X : aliased Concrete_Real_Vector := Make((2.0, 3.0, 4.0));
      X2 : constant Real_Vector_View :=
        Make(X'Access, Start => 2, Stride => 1, Length => 2);
   begin
      Assert(abs(nrm2(X)-5.3851648071345037)<Soft_Epsilon, "Nrm2 incorrect.");
      Assert(abs(nrm2(X2)-5.0)<Soft_Epsilon, "Nrm2 incorrect applied to view.");
   end Check_Nrm2;

   ----------------
   -- Check_Asum --
   ----------------

   procedure Check_Asum (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : constant Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : constant Concrete_Real_Vector := Make((-4.0, 5.0, -6.0));
   begin
      Assert(asum(X) = 6.0, "asum(X) not working");
      Assert(asum(Y) = 15.0, "asum(Y) not working (not absolute values?)");
   end Check_Asum;

   -----------------
   -- Check_Iamax --
   -----------------

   procedure Check_Iamax (T : in out Test_Cases.Test_Case'Class) is
      pragma Unreferenced (T);
      X : constant Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((4.0, -5.0, -6.0));
      Y2 : constant Real_Vector_View :=
        Make(Y'Access, Start => 1, Stride => 1, Length => 2);
   begin
      Assert(iamax(X) = 3, "iamax(X) not working");
      Assert(iamax(Y) = 3, "iamax(Y) not working");
      Assert(iamax(Y2) = 2, "iamax(Y2) not working on view");
   end Check_Iamax;

end BLAS_Real_Level1;
