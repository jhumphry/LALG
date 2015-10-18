-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 2 BLAS routines

with AUnit.Assertions; use AUnit.Assertions;

package body aBLAS_Real_Level2 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T : in out Level2_Test) is
      use AUnit.Test_Cases.Registration;
   begin

      for I of Test_Details_List loop
         Register_Routine(T, I.T, To_String(I.D));
      end loop;

   end Register_Tests;

   ------------
   -- Set_Up --
   ------------

   procedure Set_Up (T : in out Level2_Test) is
   begin
      null;
   end Set_Up;


   ----------------
   -- Check_Gemv --
   ----------------

   procedure Check_Gemv (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make(((1.0, 2.0, 3.0),
                                        (5.0, 6.0, 7.0)));
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((-4.0, 5.0));
   begin
      gemv(A, X, Y, 1.0, 2.0);
      assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X changed by GEMV operation");
      assert(A = Real_2D_Array'((1.0, 2.0, 3.0),
             (5.0, 6.0, 7.0)), "A changed by GEMV operation");
      assert(Y = Real_1D_Array'(6.0, 48.0), "Y not set correctly by GEMV operation");
      assert(gemv(A, X, 1.0) =  Real_1D_Array'(14.0, 38.0), "Function version of GEMV not working");
   end Check_Gemv;

   ----------------
   -- Check_Symv --
   ----------------

   procedure Check_Symv (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make(((1.0, 2.0),
                                               (5.0, 6.0)));
      X : aliased Concrete_Real_Vector := Make((1.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((-4.0, 5.0));
   begin
      symv(A, Upper, X, Y, 1.0, 2.0);
      assert(X = Real_1D_Array'(1.0, 3.0), "X changed by SYMV operation (upper triangular)");
      assert(A = Real_2D_Array'(
             (1.0, 2.0),
             (5.0, 6.0)), "A changed by SYMV operation (upper triangular)");
      assert(Y = Real_1D_Array'(-1.0, 30.0),
             "Y not set correctly by SYMV operation (upper triangular)");
      assert(symv(A, Upper, X, 1.0) =  Real_1D_Array'(7.0, 20.0),
             "Function version of SYMV not working (upper triangular)");

      Y := Make((-4.0, 5.0));
      symv(A, Lower, X, Y, 1.0, 2.0);
      assert(X = Real_1D_Array'(1.0, 3.0), "X changed by SYMV operation (lower triangular)");
      assert(A = Real_2D_Array'(
             (1.0, 2.0),
             (5.0, 6.0)), "A changed by SYMV operation (lower triangular)");
      assert(Y = Real_1D_Array'(8.0, 33.0),
             "Y not set correctly by SYMV operation (lower triangular)");
      assert(symv(A, Lower, X, 1.0) =  Real_1D_Array'(16.0, 23.0),
             "Function version of SYMV not working (lower triangular)");

   end Check_Symv;

   ---------------
   -- Check_Ger --
   ---------------

   procedure Check_Ger (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((-4.0, 5.0));
      A : aliased Concrete_Real_Matrix := Zeros(3,2);
   begin
      ger(X, Y, A, 2.0);
      assert(X = Real_1D_Array'(1.0, 2.0, 3.0), "X changed by GER operation");
      assert(Y = Real_1D_Array'(-4.0, 5.0), "Y changed by GER operation");
      assert(A = Real_2D_Array'(((-8.0, 10.0),
             (-16.0, 20.0),
             (-24.0, 30.0))), "A not set correctly by GER operation");

   end Check_Ger;

   --------------------
   -- Check_Syr_Syr2 --
   --------------------

   procedure Check_Syr_Syr2 (T : in out Test_Cases.Test_Case'Class) is
      X : aliased Concrete_Real_Vector := Make((1.0, 2.0, 3.0));
      Y : aliased Concrete_Real_Vector := Make((-1.0, 2.0, -1.0));
      A : aliased Concrete_Real_Matrix := Zeros(3,3);
   begin
      syr(X, A, Upper, 1.0);
      assert(X = Real_1D_Array'(1.0, 2.0, 3.0),
             "X changed by SYR operation (upper triangular)");
      assert(A = Real_2D_Array'((
             (1.0, 2.0, 3.0),
             (0.0, 4.0, 6.0),
             (0.0, 0.0, 9.0))),
             "A not set correctly by SYR operation (upper triangular)");

      syr2(X, Y, A, Lower, 2.0);
      assert(X = Real_1D_Array'(1.0, 2.0, 3.0),
             "X changed by SYR operation (lower triangular)");
      assert(Y = Real_1D_Array'(-1.0, 2.0, -1.0),
             "Y changed by SYR operation (lower triangular)");
      assert(A = Real_2D_Array'((
             (-3.0, 2.0, 3.0),
             ( 0.0,20.0, 6.0),
             (-8.0, 8.0,-3.0))),
             "A not set correctly by SYR operation (upper triangular)");

   end Check_Syr_Syr2;

end aBLAS_Real_Level2;
