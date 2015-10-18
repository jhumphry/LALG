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

    ----------------
   -- Check_Symm --
   ----------------

   procedure Check_Symm (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make((
                                               (1.0, 3.0),
                                               (5.0, 7.0)
                                              ));
      B : aliased Concrete_Real_Matrix := Make((
                                               (11.0, 13.0),
                                               (17.0, 19.0)
                                              ));
      C : aliased Concrete_Real_Matrix := Identity(2);

   begin

      -- Left, Upper
      symm(A     => A,
           SIDE  => Left,
           UPLO  => Upper,
           B     => B,
           C     => C,
           ALPHA => 1.0,
           BETA  => 2.0);

      assert(A = Real_2D_Array'(
             (1.0, 3.0),
             (5.0, 7.0)
            ), "A changed by SYMM operation (Left, Upper)");

      assert(B = Real_2D_Array'(
             (11.0, 13.0),
             (17.0, 19.0)
             ), "B changed by SYMM operation (Left, Upper)");

      assert(C = Real_2D_Array'(
             ( 64.0,  70.0),
             (152.0, 174.0)
             ), "C not set correctly by SYMM operation (Left, Upper)");

      -- Right, Upper
      C := Identity(2);

      symm(A     => A,
           SIDE  => Right,
           UPLO  => Upper,
           B     => B,
           C     => C,
           ALPHA => 1.0,
           BETA  => 2.0);

      assert(C = Real_2D_Array'(
             ( 52.0, 124.0),
             ( 74.0, 186.0)
            ), "C not set correctly by SYMM operation (Right, Upper)");

      -- Left, Lower
      C := Identity(2);

      symm(A     => A,
           SIDE  => Left,
           UPLO  => Lower,
           B     => B,
           C     => C,
           ALPHA => 1.0,
           BETA  => 2.0);

      assert(C = Real_2D_Array'(
             ( 98.0, 108.0),
             (174.0, 200.0)
            ), "C not set correctly by SYMM operation (Left, Lower)");

      -- Right, Lower
      C := Identity(2);

      symm(A     => A,
           SIDE  => Right,
           UPLO  => Lower,
           B     => B,
           C     => C,
           ALPHA => 1.0,
           BETA  => 2.0);

      assert(C = Real_2D_Array'(
             ( 78.0, 146.0),
             (112.0, 220.0)
            ), "C not set correctly by SYMM operation (Left, Lower)");

      assert(symm(A     => A,
                  SIDE  => Left,
                  UPLO  => Upper,
                  B     => B,
                  ALPHA => 1.0) =
               Real_2D_Array'(
                 ( 62.0, 70.0),
                 (152.0, 172.0)
                ), "Function version of SYMM not working");
   end Check_Symm;

   ----------------
   -- Check_Syrk --
   ----------------

   procedure Check_Syrk (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make((
                                               (1.0, 3.0),
                                               (5.0, 7.0)
                                              ));
      C_Original : constant Concrete_Real_Matrix := Make((
                                               (11.0, 13.0),
                                               (17.0, 19.0)
                                              ));
      C : aliased Concrete_Real_Matrix := C_Original;

   begin

      -- No, Upper
      syrk(A     => A,
           TRANS => No_Transpose,
           UPLO  => Upper,
           C     => C,
           ALPHA => 2.0,
           BETA  => 1.0);

      assert(A = Real_2D_Array'(
             (1.0, 3.0),
             (5.0, 7.0)
            ), "A changed by SYRK operation (No, Upper)");

      assert(C = Real_2D_Array'(
             (31.0,  65.0),
             (17.0, 167.0)
             ), "C not set correctly by SYRK operation (No, Upper)");

      -- Transpose, Upper
      C := C_Original;

      syrk(A     => A,
           TRANS => Transpose,
           UPLO  => Upper,
           C     => C,
           ALPHA => 2.0,
           BETA  => 1.0);

      assert(C = Real_2D_Array'(
             (63.0,  89.0),
             (17.0, 135.0)
             ), "C not set correctly by SYRK operation (Transpose, Upper)");

      -- No, Lower
      C := C_Original;

      syrk(A     => A,
           TRANS => No_Transpose,
           UPLO  => Lower,
           C     => C,
           ALPHA => 2.0,
           BETA  => 1.0);

      assert(C = Real_2D_Array'(
             (31.0,  13.0),
             (69.0, 167.0)
             ), "C not set correctly by SYRK operation (No, Lower)");

      -- Transpose, Lower
      C := C_Original;

      syrk(A     => A,
           TRANS => Transpose,
           UPLO  => Lower,
           C     => C,
           ALPHA => 2.0,
           BETA  => 1.0);

      assert(C = Real_2D_Array'(
             (63.0,  13.0),
             (93.0, 135.0)
             ), "C not set correctly by SYRK operation (Transpose, Lower)");

      assert(syrk(A     => A,
                  TRANS => No_Transpose,
                  UPLO  => Upper,
                  ALPHA => 2.0) =
               Real_2D_Array'(
                 ( 20.0, 52.0),
                 (  0.0, 148.0)
                ), "Function version of SYRK not working");
   end Check_Syrk;

end aBLAS_Test_Suite.Real_Level3;
