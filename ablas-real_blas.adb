-- aBLAS
-- An Ada 2012 binding to BLAS

with aBLAS.Real_BLAS.Imports;

with Interfaces.Fortran;

package body aBLAS.Real_BLAS is

   package Fortran_Imports is new aBLAS.Real_BLAS.Imports;
   use Fortran_Imports;

   function TF(Item : in Character) return Interfaces.Fortran.Character_Set
               renames Interfaces.Fortran.To_Fortran;

   type Map_Trans_Op_Array is array (Trans_Op) of IntFort.Character_Set;
   Map_Trans_Op : constant Map_Trans_Op_Array := (No_Transpose => TF('N'),
                                                  Transpose => TF('T'),
                                                  Conj_Transpose => TF('C'));

   type Map_UpLo_Part_Array is array (UpLo_Part) of IntFort.Character_Set;
   Map_UpLo_Part : constant Map_UpLo_Part_Array := (Upper => TF('U'),
                                                    Lower => TF('L'));

   type Map_Diag_Unit_Array is array (Diag_Unit) of IntFort.Character_Set;
   Map_Diag_Unit : constant Map_Diag_Unit_Array := (Non_Unit_Diag => TF('N'),
                                                    Unit_Diag => TF('D'));

   type Map_Side_Op_Array is array (Side_Op) of IntFort.Character_Set;
   Map_Side_Op : constant Map_Side_Op_Array := (Left => TF('L'),
                                                Right => TF('R'));

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

   ---------
   -- syr --
   ---------

   procedure syr(X : in Real_Vector'Class;
                 A : in out Real_Matrix'Class;
                 UPLO : in UpLo_Part;
                 ALPHA : in Real := 1.0)
   is
   begin
      case Precision is
         when Single =>
            SSYR(UPLO => Map_UpLo_Part(UPLO),
                 N => FP(A.Columns),
                 ALPHA => ALPHA,
                 X => X.Constant_Handle,
                 INCX => FP(X.Stride),
                 A => A.Handle,
                 LDA => FP(A.Leading_Dimension));
         when Double =>
            DSYR(UPLO => Map_UpLo_Part(UPLO),
                 N => FP(A.Columns),
                 ALPHA => ALPHA,
                 X => X.Constant_Handle,
                 INCX => FP(X.Stride),
                 A => A.Handle,
                 LDA => FP(A.Leading_Dimension));
      end case;
   end syr;

   ----------
   -- syr2 --
   ----------

   procedure syr2(X : in Real_Vector'Class;
                  Y : in Real_Vector'Class;
                  A : in out Real_Matrix'Class;
                  UPLO : in UpLo_Part;
                  ALPHA : in Real := 1.0)
   is
   begin
      case Precision is
         when Single =>
            SSYR2(UPLO => Map_UpLo_Part(UPLO),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  Y => Y.Constant_Handle,
                  INCY => FP(Y.Stride),
                  A => A.Handle,
                  LDA => FP(A.Leading_Dimension));
         when Double =>
            DSYR2(UPLO => Map_UpLo_Part(UPLO),
                  N => FP(A.Columns),
                  ALPHA => ALPHA,
                  X => X.Constant_Handle,
                  INCX => FP(X.Stride),
                  Y => Y.Constant_Handle,
                  INCY => FP(Y.Stride),
                  A => A.Handle,
                  LDA => FP(A.Leading_Dimension));
      end case;
   end syr2;

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

   procedure symm(A : in Real_Matrix'Class;
                  SIDE : in Side_Op;
                  UPLO : in UpLo_Part;
                  B : in Real_Matrix'Class;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0) is
   begin
      case Precision is
         when Single =>
            SSYMM(SIDE => Map_Side_Op(SIDE),
                  UPLO => Map_UpLo_Part(UPLO),
                  M => FP(C.Rows),
                  N => FP(C.Columns),
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
            DSYMM(SIDE => Map_Side_Op(SIDE),
                  UPLO => Map_UpLo_Part(UPLO),
                  M => FP(C.Rows),
                  N => FP(C.Columns),
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
   end symm;

   function symm(A : in Real_Matrix'Class;
                 SIDE : in Side_Op;
                 UPLO : in UpLo_Part;
                 B : in Real_Matrix'Class;
                 ALPHA : in Real := 1.0)
                 return Concrete_Real_Matrix is
      C : Concrete_Real_Matrix(M => (case SIDE is
                                        when Left => A.Rows,
                                        when Right => B.Rows),
                               N => (case SIDE is
                                        when Left => B.Columns,
                                        when Right => A.Columns));
      -- As Beta is being set to zero, C should only be written to and not read
      -- so it does not matter that it is uninitialised.
   begin
      symm(A     => A,
           SIDE  => SIDE,
           UPLO  => UPLO,
           B     => B,
           C     => C,
           ALPHA => ALPHA,
           BETA  => 0.0);
      return C;
   end symm;

   procedure syrk(A : in Real_Matrix'Class;
                  TRANS : in Real_Trans_Op;
                  UPLO : in UpLo_Part;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0) is
      K : FP := (if TRANS = No_Transpose then
                    FP(A.Columns)
                 else
                    FP(A.Rows));
   begin
      case Precision is
         when Single =>
            SSYRK(UPLO => Map_UpLo_Part(UPLO),
                  TRANS => Map_Trans_Op(TRANS),
                  N => FP(C.Rows),
                  K => K,
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  BETA => BETA,
                  C => C.Handle,
                  LDC => FP(C.Leading_Dimension)
                 );
         when Double =>
            DSYRK(UPLO => Map_UpLo_Part(UPLO),
                  TRANS => Map_Trans_Op(TRANS),
                  N => FP(C.Rows),
                  K => K,
                  ALPHA => ALPHA,
                  A => A.Constant_Handle,
                  LDA => FP(A.Leading_Dimension),
                  BETA => BETA,
                  C => C.Handle,
                  LDC => FP(C.Leading_Dimension)
                 );
      end case;
   end syrk;

   function syrk(A : in Real_Matrix'Class;
                 TRANS : in Real_Trans_Op;
                 UPLO : in UpLo_Part;
                 ALPHA : in Real := 1.0)
                 return Concrete_Real_Matrix is
      C_Order : Positive := (case TRANS is
                                when No_Transpose => A.Rows,
                                when Transpose => A.Columns);
      C : Concrete_Real_Matrix := Zeros(Rows => C_Order,
                                        Columns => C_Order);
   begin
      syrk(A     => A,
           TRANS => TRANS,
           UPLO  => UPLO,
           C     => C,
           ALPHA => ALPHA,
           BETA  => 0.0);
      return C;
   end syrk;

   procedure syr2k(A : in Real_Matrix'Class;
                   B : in Real_Matrix'Class;
                   TRANS : in Real_Trans_Op;
                   UPLO : in UpLo_Part;
                   C : in out Real_Matrix'Class;
                   ALPHA : in Real := 1.0;
                   BETA : in Real := 0.0) is
      K : FP := (if TRANS = No_Transpose then
                    FP(A.Columns)
                 else
                    FP(A.Rows));
   begin
      case Precision is
         when Single =>
            SSYR2K(UPLO => Map_UpLo_Part(UPLO),
                   TRANS => Map_Trans_Op(TRANS),
                   N => FP(C.Rows),
                   K => K,
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
            DSYR2K(UPLO => Map_UpLo_Part(UPLO),
                   TRANS => Map_Trans_Op(TRANS),
                   N => FP(C.Rows),
                   K => K,
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
   end syr2k;

   function syr2k(A : in Real_Matrix'Class;
                  B : in Real_Matrix'Class;
                  TRANS : in Real_Trans_Op;
                  UPLO : in UpLo_Part;
                  ALPHA : in Real := 1.0)
                  return Concrete_Real_Matrix is
      C_Order : Positive := (case TRANS is
                                when No_Transpose => A.Rows,
                                when Transpose => A.Columns);
      C : Concrete_Real_Matrix := Zeros(Rows => C_Order,
                                        Columns => C_Order);
   begin
      syr2k(A    => A,
            B    => B,
            TRANS => TRANS,
            UPLO  => UPLO,
            C     => C,
            ALPHA => ALPHA,
            BETA  => 0.0);
      return C;
   end syr2k;

end aBLAS.Real_BLAS;
