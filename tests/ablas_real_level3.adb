-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests for Level 3 BLAS routines

with AUnit.Assertions; use AUnit.Assertions;

package body aBLAS_Real_Level3 is

   --------------------
   -- Register_Tests --
   --------------------

   procedure Register_Tests (T : in out Level3_Test) is
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

   begin
      gemm(A => A,
           B => B,
           C => C,
           ALPHA => 1.0,
           BETA => 2.0,
           TRANA => No_Transpose,
           TRANB => No_Transpose);
      Assert(A = Real_2D_Array'(
             (1.0, 2.0),
             (3.0, 4.0),
             (5.0, 6.0)
            ), "A changed by GEMM operation");

      Assert(B = Real_2D_Array'(
             (1.0, 2.0, 1.0),
             (2.0, 1.0, 2.0)), "B changed by GEMM operation");

      Assert(C = Real_2D_Array'(
             ( 7.0,  4.0,  5.0),
             (11.0, 12.0, 11.0),
             (17.0, 16.0, 19.0)), "C not set correctly by GEMM operation");

      Assert(gemm(A, B, 1.0) = Real_2D_Array'(
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

      Assert(A = Real_2D_Array'(
             (1.0, 3.0),
             (5.0, 7.0)
            ), "A changed by SYMM operation (Left, Upper)");

      Assert(B = Real_2D_Array'(
             (11.0, 13.0),
             (17.0, 19.0)
            ), "B changed by SYMM operation (Left, Upper)");

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
             ( 78.0, 146.0),
             (112.0, 220.0)
            ), "C not set correctly by SYMM operation (Left, Lower)");

      Assert(symm(A     => A,
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

      Assert(A = Real_2D_Array'(
             (1.0, 3.0),
             (5.0, 7.0)
            ), "A changed by SYRK operation (No, Upper)");

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
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

      Assert(C = Real_2D_Array'(
             (63.0,  13.0),
             (93.0, 135.0)
            ), "C not set correctly by SYRK operation (Transpose, Lower)");

      Assert(syrk(A     => A,
                  TRANS => No_Transpose,
                  UPLO  => Upper,
                  ALPHA => 2.0) =
               Real_2D_Array'(
                 ( 20.0, 52.0),
                 (  0.0, 148.0)
                ), "Function version of SYRK not working");
   end Check_Syrk;

   -----------------
   -- Check_Syr2k --
   -----------------

   procedure Check_Syr2k (T : in out Test_Cases.Test_Case'Class) is
      A : aliased Concrete_Real_Matrix := Make((
                                               (1.0, 3.0),
                                               (5.0, 7.0)
                                              ));
      B : aliased Concrete_Real_Matrix := Make((
                                               (11.0, 13.0),
                                               (17.0, 19.0)
                                              ));
      C_Original : constant Concrete_Real_Matrix := Make((
                                                         (23.0, 29.0),
                                                         (31.0, 37.0)
                                                        ));
      C : aliased Concrete_Real_Matrix := C_Original;

   begin

      -- No, Upper
      syr2k(A     => A,
            B     => B,
            TRANS => No_Transpose,
            UPLO  => Upper,
            C     => C,
            ALPHA => 2.0,
            BETA  => 1.0);

      Assert(A = Real_2D_Array'(
             (1.0, 3.0),
             (5.0, 7.0)
            ), "A changed by SYR2K operation (No, Upper)");

      Assert(B = Real_2D_Array'(
             (11.0, 13.0),
             (17.0, 19.0)
            ), "B changed by SYR2K operation (No, Upper)");

      Assert(C = Real_2D_Array'(
             (223.0, 469.0),
             ( 31.0, 909.0)
            ), "C not set correctly by SYR2K operation (No, Upper)");

      -- Transpose, Upper
      C := C_Original;

      syr2k(A     => A,
            B     => B,
            TRANS => Transpose,
            UPLO  => Upper,
            C     => C,
            ALPHA => 2.0,
            BETA  => 1.0);

      Assert(C = Real_2D_Array'(
             (407.0, 549.0),
             ( 31.0, 725.0)
            ), "C not set correctly by SYR2K operation (Transpose, Upper)");

      -- No, Lower
      C := C_Original;

      syr2k(A     => A,
            B     => B,
            TRANS => No_Transpose,
            UPLO  => Lower,
            C     => C,
            ALPHA => 2.0,
            BETA  => 1.0);

      Assert(C = Real_2D_Array'(
             (223.0,  29.0),
             (471.0, 909.0)
            ), "C not set correctly by SYR2K operation (No, Lower)");

      -- Transpose, Lower
      C := C_Original;

      syr2k(A     => A,
            B     => B,
            TRANS => Transpose,
            UPLO  => Lower,
            C     => C,
            ALPHA => 2.0,
            BETA  => 1.0);

      Assert(C = Real_2D_Array'(
             (407.0,  29.0),
             (551.0, 725.0)
            ), "C not set correctly by SYR2K operation (Transpose, Lower)");

      Assert(syr2k(A     => A,
                   B     => B,
                   TRANS => No_Transpose,
                   UPLO  => Upper,
                   ALPHA => 2.0) =
               Real_2D_Array'(
                 (200.0, 440.0),
                 (  0.0, 872.0)
                ), "Function version of SYR2K not working");
   end Check_Syr2k;

end aBLAS_Real_Level3;