-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 2 BLAS routines

with AUnit.Assertions; use AUnit.Assertions;

with aBLAS_Test_Suite;

package body aBLAS_Test_Suite.Real_Level2 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T: in out Level2_Test) is
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

end aBLAS_Test_Suite.Real_Level2;
