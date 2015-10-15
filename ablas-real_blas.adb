-- aBLAS
-- An Ada 2012 binding to BLAS

with aBLAS.Real_BLAS.Imports;

package body aBLAS.Real_BLAS is

   package Fortran_Imports is new aBLAS.Real_BLAS.Imports;
   use Fortran_Imports;

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   ----------
   -- rotg --
   ----------

   procedure rotg(a, b : in out Real; c, s : out Real) is
   begin
       case Precision is
         when Single => SROTG(A,B,C,S);
         when Double => DROTG(A,B,C,S);
      end case;
   end rotg;

   -----------
   -- rotmg --
   -----------

   procedure rotmg(d1, d2 : in out Real;
                   x1 : in out Real;
                   y1 : in Real;
                   params : out Modified_Givens_Params) is
   begin
       case Precision is
         when Single => SROTMG(D1, D2, X1, Y1, PARAMS);
         when Double => DROTMG(D1, D2, X1, Y1, PARAMS);
      end case;
   end rotmg;

   ---------
   -- rot --
   ---------

   procedure rot(X : in out Real_Vector'Class;
                 Y : in out Real_Vector'Class;
                 C : in Real;
                 S : in Real)
   is
   begin
      case Precision is
         when Single => SROT(N => FP(X.Length),
                             SX => X.Handle,
                             INCX => FP(X.Stride),
                             SY => Y.Handle,
                             INCY => FP(Y.Stride),
                             C => C,
                             S => S);
         when Double => DROT(N => FP(X.Length),
                             DX => X.Handle,
                             INCX => FP(X.Stride),
                             DY => Y.Handle,
                             INCY => FP(Y.Stride),
                             C => C,
                             S => S);
      end case;
   end rot;

   ----------
   -- rotm --
   ----------

   procedure rotm(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  PARAMS : in out Modified_Givens_Params)
   is
   begin
      case Precision is
         when Single => SROTM(N => FP(X.Length),
                              SX => X.Handle,
                              INCX => FP(X.Stride),
                              SY => Y.Handle,
                              INCY => FP(Y.Stride),
                              PARAMS => PARAMS);
         when Double => DROTM(N => FP(X.Length),
                              DX => X.Handle,
                              INCX => FP(X.Stride),
                              DY => Y.Handle,
                              INCY => FP(Y.Stride),
                              PARAMS => PARAMS);
      end case;
   end rotm;

   ----------
   -- swap --
   ----------

   procedure swap(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class)
    is
   begin
      case Precision is
         when Single => SSWAP(N => FP(X.Length),
                              SX => X.Handle,
                              INCX => FP(X.Stride),
                              SY => Y.Handle,
                              INCY => FP(Y.Stride));
         when Double => DSWAP(N => FP(X.Length),
                              DX => X.Handle,
                              INCX => FP(X.Stride),
                              DY => Y.Handle,
                              INCY => FP(Y.Stride));
      end case;
   end swap;

   ----------
   -- scal --
   ----------

   procedure scal(X : in out Real_Vector'Class;
                  A : in Real := 1.0)
    is
   begin
      case Precision is
         when Single => SSCAL(N => FP(X.Length),
                              SA => A,
                              SX => X.Handle,
                              INCX =>  FP(X.Stride));
         when Double => DSCAL(N =>  FP(X.Length),
                              DA => A,
                              DX => X.Handle,
                              INCX =>  FP(X.Stride));
      end case;
   end scal;

   ----------
   -- copy --
   ----------

   procedure copy(X : in Real_Vector'Class;
                  Y : out Real_Vector'Class)
    is
   begin
      case Precision is
         when Single => SCOPY(N => FP(X.Length),
                              SX => X.Constant_Handle,
                              INCX => FP(X.Stride),
                              SY => Y.Handle,
                              INCY => FP(Y.Stride));
         when Double => DCOPY(N => FP(X.Length),
                              DX => X.Constant_Handle,
                              INCX => FP(X.Stride),
                              DY => Y.Handle,
                              INCY => FP(Y.Stride));
      end case;
   end copy;

   ----------
   -- axpy --
   ----------

   procedure axpy(X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  A : in Real := 1.0)
   is
   begin
      case Precision is
         when Single => SAXPY(N => FP(X.Length),
                              SA => A,
                              SX => X.Constant_Handle,
                              INCX => FP(X.Stride),
                              SY => Y.Handle,
                              INCY => FP(Y.Stride));
         when Double => DAXPY(N => FP(X.Length),
                              DA => A,
                              DX => X.Constant_Handle,
                              INCX => FP(X.Stride),
                              DY => Y.Handle,
                              INCY => FP(Y.Stride));
      end case;
   end axpy;

   ---------
   -- dot --
   ---------

   function dot
     (X, Y : in Real_Vector'Class)
      return Real
   is
   begin
      case Precision is
         when Single => return SDOT(N => FP(X.Length),
                                    SX => X.Constant_Handle,
                                    INCX => FP(X.Stride),
                                    SY => Y.Constant_Handle,
                                    INCY => FP(Y.Stride));
         when Double => return DDOT(N => FP(X.Length),
                                    DX => X.Constant_Handle,
                                    INCX => FP(X.Stride),
                                    DY => Y.Constant_Handle,
                                    INCY => FP(Y.Stride));
      end case;
   end dot;

   ------------
   -- sdsdot --
   ------------

   function sdsdot
     (X, Y: in Real_Vector'Class;
      B : in Real := 0.0)
      return Real
   is
   begin
      case Precision is
         when Single => return SDSDOT(N => FP(X.Length),
                                      SB => B,
                                      SX => X.Constant_Handle,
                                      INCX => FP(X.Stride),
                                      SY => Y.Constant_Handle,
                                      INCY => FP(Y.Stride));
         when Double => return DDOT(N => FP(X.Length),
                                    DX => X.Constant_Handle,
                                    INCX => FP(X.Stride),
                                    DY => Y.Constant_Handle,
                                    INCY => FP(Y.Stride)) + B;
      end case;
   end sdsdot;

   ----------
   -- nrm2 --
   ----------

   function nrm2
     (X : in Real_Vector'Class)
      return Real
   is
   begin
      case Precision is
         when Single => return SNRM2(N => FP(X.Length),
                                     X => X.Constant_Handle,
                                     INCX => FP(X.Stride));
         when Double => return DNRM2(N => FP(X.Length),
                                     X => X.Constant_Handle,
                                     INCX => FP(X.Stride));
      end case;
   end nrm2;

   ----------
   -- asum --
   ----------

   function asum(X : in Real_Vector'Class) return Real is
   begin
      case Precision is
         when Single =>
            return SASUM(N => FP(X.Length),
                         SX => X.Constant_Handle,
                         INCX => FP(X.Stride));
         when Double=>
            return DASUM(N => FP(X.Length),
                         DX => X.Constant_Handle,
                         INCX => FP(X.Stride));
      end case;
   end asum;

   function iamax(X : in Real_Vector'Class) return Integer
   is
   begin
      case Precision is
         when Single => return Natural(ISAMAX(N => FP(X.Length),
                                              SX => X.Constant_Handle,
                                              INCX => FP(X.Stride)));
         when Double => return Natural(IDAMAX(N => FP(X.Length),
                                              DX => X.Constant_Handle,
                                              INCX => FP(X.Stride)));
      end case;
   end iamax;

   -- *************
   -- *************
   -- ** Level 2 **
   -- *************
   -- *************

   ----------
   -- gemv --
   ----------

   procedure gemv(A : in Real_Matrix'Class;
                  X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANS : in Real_Trans_Op := No_Transpose) is
   begin
      case Precision is
         when Single =>
            SGEMV(TRANS => Map_Trans_Op(TRANS),
                  M => FP(A.Rows),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  BETA => BETA,
                  Y => Y.Handle,
                  INCY => FP(Y.Stride));
         when Double =>
            DGEMV(TRANS => Map_Trans_Op(TRANS),
                  M => FP(A.Rows),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  BETA => BETA,
                  Y => Y.Handle,
                  INCY => FP(Y.Stride));
      end case;
   end gemv;

   function gemv(A : in Real_Matrix'Class;
                 X : in Real_Vector'Class;
                 ALPHA : in Real := 1.0;
                 TRANS : in Real_Trans_Op := No_Transpose)
                 return Real_Vector'Class is
      Y : Concrete_Real_Vector(N => A.Rows);
      -- As Beta is being set to zero, Y should only be written to and not read
      -- so it does not matter that it is uninitialised.
   begin
      gemv(A => A,
           X => X,
           Y => Y,
           ALPHA => ALPHA,
           BETA => 0.0,
           TRANS => TRANS);
      return Y;
   end gemv;

   ----------
   -- symv --
   ----------

   procedure symv(A : in Real_Matrix'Class;
                  UPLO : in UpLo_Part;
                  X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0) is
   begin
      case Precision is
         when Single =>
            SSYMV(UPLO => Map_UpLo_Part(UPLO),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  BETA => BETA,
                  Y => Y.Handle,
                  INCY => FP(Y.Stride));
         when Double =>
            DSYMV(UPLO => Map_UpLo_Part(UPLO),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  BETA => BETA,
                  Y => Y.Handle,
                  INCY => FP(Y.Stride));
      end case;
   end symv;

   function symv(A : in Real_Matrix'Class;
                 UPLO : in UpLo_Part;
                 X : in Real_Vector'Class;
                 ALPHA : in Real := 1.0)
                 return Real_Vector'Class is
      Y : Concrete_Real_Vector(N => A.Rows);
      -- As Beta is being set to zero, Y should only be written to and not read
      -- so it does not matter that it is uninitialised.
   begin
      symv(A => A,
           UPLO => UPLO,
           X => X,
           Y => Y,
           ALPHA => ALPHA,
           BETA => 0.0);
      return Y;
   end symv;

   ---------
   -- ger --
   ---------

   procedure ger(X : in Real_Vector'Class;
                 Y : in Real_Vector'Class;
                 A : in out Real_Matrix'Class;
                 ALPHA : in Real := 1.0)
   is
   begin
      case Precision is
         when Single =>
            SGER(M => FP(A.Rows),
                 N => FP(A.Columns),
                 ALPHA => ALPHA,
                 X => X.Constant_Handle,
                 INCX => FP(X.Stride),
                 Y => Y.Constant_Handle,
                 INCY => FP(Y.Stride),
                 A => A.Handle,
                 LDA => FP(A.Leading_Dimension));
         when Double =>
            DGER(M => FP(A.Rows),
                 N => FP(A.Columns),
                 ALPHA => ALPHA,
                 X => X.Constant_Handle,
                 INCX => FP(X.Stride),
                 Y => Y.Constant_Handle,
                 INCY => FP(Y.Stride),
                 A => A.Handle,
                 LDA => FP(A.Leading_Dimension));
      end case;
   end ger;

   -- *************
   -- *************
   -- ** Level 3 **
   -- *************
   -- *************

   procedure gemm(A : in Real_Matrix'Class;
                  B : in Real_Matrix'Class;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANA : in Real_Trans_Op := No_Transpose;
                  TRANB : in Real_Trans_Op := No_Transpose) is
   begin
      case Precision is
         when Single =>
            SGEMM(TRANA => Map_Trans_Op(TRANA),
                  TRANB => Map_Trans_Op(TRANB),
                  M => FP(A.Rows),
                  N => FP(B.Columns),
                  K => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  B => B.Constant_Handle,
                  LDB => FP(B.Leading_Dimension),
                  BETA => BETA,
                  C => C.Handle,
                  LDC => FP(C.Leading_Dimension)
                 );
         when Double =>
            DGEMM(TRANA => Map_Trans_Op(TRANA),
                  TRANB => Map_Trans_Op(TRANB),
                  M => FP(A.Rows),
                  N => FP(B.Columns),
                  K => FP(A.Columns),
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  B => B.Constant_Handle,
                  LDB => FP(B.Leading_Dimension),
                  BETA => BETA,
                  C => C.Handle,
                  LDC => FP(C.Leading_Dimension)
                 );
      end case;
   end gemm;

   function gemm(A : in Real_Matrix'Class;
                 B : in Real_Matrix'Class;
                 ALPHA : in Real := 1.0;
                 TRANA : in Real_Trans_Op := No_Transpose;
                 TRANB : in Real_Trans_Op := No_Transpose)
                 return Concrete_Real_Matrix is
      C : Concrete_Real_Matrix(A.Rows, B.Columns);
      -- As Beta is being set to zero, C should only be written to and not read
      -- so it does not matter that it is uninitialised.
   begin
      gemm(A => A,
           B => B,
           C => C,
           ALPHA => ALPHA,
           BETA => 0.0,
           TRANA => TRANA,
           TRANB => TRANB);
      return C;
   end gemm;

end aBLAS.Real_BLAS;
