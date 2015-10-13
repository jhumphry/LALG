-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 3 BLAS routines

with AUnit.Assertions; use AUnit.Assertions;

with aBLAS_Test_Suite;

package body aBLAS_Test_Suite.Real_Level3 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T: in out Level3_Test) is
      use AUnit.Test_Cases.Registration;
   begin

      for I of Test_Details_List loop
         Register_Routine(T, I.T, To_String(I.D));
      end loop;

   end Register_Tests;

   ------------
   -- Set_Up --
   ------------

   procedure Set_Up (T : in out Level3_Test) is
   begin
      null;
   end Set_Up;

   ----------------
   -- Check_Gemm --
   ----------------

   procedure Check_Gemm (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make((
                                               (1.0, 2.0),
                                               (3.0, 4.0),
                                               (5.0, 6.0)
                                              ));
      B : aliased Concrete_Real_Matrix := Make((
                                               (1.0, 2.0, 1.0),
                                               (2.0, 1.0, 2.0)
                                              ));
      C : aliased Concrete_Real_Matrix := Identity(3);
      Expected_C : aliased Concrete_Real_Matrix := Make((
                                               ( 7.0,  4.0,  5.0),
                                               (11.0, 12.0, 11.0),
                                               (17.0, 16.0, 19.0)
                                              ));
   begin
      gemm(A => A,
           B => B,
           C => C,
           ALPHA => 1.0,
           BETA => 2.0,
           TRANA => No_Transpose,
           TRANB => No_Transpose);
      assert(A = Real_2D_Array'(
             (1.0, 2.0),
             (3.0, 4.0),
             (5.0, 6.0)
            ), "A changed by GEMM operation");

      assert(B = Real_2D_Array'(
             (1.0, 2.0, 1.0),
             (2.0, 1.0, 2.0)), "B changed by GEMM operation");

      assert(C = Real_2D_Array'(
             ( 7.0,  4.0,  5.0),
             (11.0, 12.0, 11.0),
             (17.0, 16.0, 19.0)), "C not set correctly by GEMM operation");

      assert(gemm(A, B, 1.0) = Real_2D_Array'(
             ( 5.0,  4.0,  5.0),
             (11.0, 10.0, 11.0),
             (17.0, 16.0, 17.0)), "Function version of GEMM not working");
   end Check_Gemm;

end aBLAS_Test_Suite.Real_Level3;
