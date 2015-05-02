-- aBLAS
-- An Ada 2012 binding to BLAS

with aBLAS.Real_BLAS.Imports;
with aBLAS.Internal;

package body aBLAS.Real_BLAS is

   package Fortran_Imports is new aBLAS.Real_BLAS.Imports;
   use Fortran_Imports;

   Map_Trans_Op : Internal.Map_Trans_Op_Array renames Internal.Map_Trans_Op;

   ---------
   -- dot --
   ---------

   function dot
     (X, Y: in Real_Vector;
      INCX : in Increment := 1;
      INCY : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then X'Length else N);
   begin
      case Precision is
         when Single => return SDOT(NN,X,INCX,Y,INCY);
         when Double => return DDOT(NN,X,INCX,Y,INCY);
      end case;
   end dot;

   ------------
   -- sdsdot --
   ------------

   function sdsdot
     (SX, SY: in Real_Vector;
      SB : Real := 0.0;
      INCX : in Increment := 1;
      INCY : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then SX'Length else N);
   begin
      case Precision is
         when Single => return SDSDOT(NN,SB,SX,INCX,SY,INCY);
         when Double => return DDOT(NN,SX,INCX,SY,INCY) + SB;
      end case;
   end sdsdot;

   ----------
   -- nrm2 --
   ----------

   function nrm2
     (X : in Real_Vector;
      INCX : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then X'Length else N);
   begin
      case Precision is
         when Single => return SNRM2(NN,X,INCX);
         when Double => return DNRM2(NN,X,INCX);
      end case;
   end nrm2;

   ----------
   -- asum --
   ----------

   function asum
     (X : in Real_Vector;
      INCX : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then X'Length else N);
   begin
      case Precision is
         when Single => return SASUM(NN,X,INCX);
         when Double => return DASUM(NN,X,INCX);
      end case;
   end asum;

   ----------
   -- gemv --
   ----------

   procedure gemv(A : in Real_Matrix;
                  X: in Real_Vector;
                  Y : in out Real_Vector;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANS : in Real_Trans_Op := No_Transpose;
                  INCX, INCY : in Increment := 1;
                  M, N : in Vector_Size := 0;
                  Convention : in Matrix_Convention := Default_Matrix_Convention) is
      TRANS_P : Real_Trans_Op;
      M_P : constant Vector_Size := (if M = 0 then A'Length(2) else M);
      N_P : constant Vector_Size := (if N = 0 then A'Length(1) else N);
   begin
      case Convention is
         when Row_Major =>
            TRANS_P := (if TRANS = No_Transpose then Transpose else No_Transpose);
         when Column_Major =>
            TRANS_P := TRANS;
      end case;

      case Precision is
         when Single =>
            SGEMV(TRANS => Map_Trans_Op(TRANS_P),
                  M => M_P,
                  N => N_P,
                  ALPHA => ALPHA,
                  A => A,
                  LDA => A'Length(2),
                  X => X,
                  INCX => INCX,
                  BETA => BETA,
                  Y => Y,
                  INCY => INCY);
         when Double =>
            DGEMV(TRANS => Map_Trans_Op(TRANS_P),
                  M => M_P,
                  N => N_P,
                  ALPHA => ALPHA,
                  A => A,
                  LDA => A'Length(2),
                  X => X,
                  INCX => INCX,
                  BETA => BETA,
                  Y => Y,
                  INCY => INCY);
      end case;
   end gemv;

   -- *************
   -- *************
   -- ** Level 3 **
   -- *************
   -- *************

   procedure gemm(A : in Real_Matrix;
                  B : in Real_Matrix;
                  C : in out Real_Matrix;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANA : in Real_Trans_Op := No_Transpose;
                  TRANB : in Real_Trans_Op := No_Transpose;
                  M, N, K : in Vector_Size := 0;
                  Convention : in Matrix_Convention := Default_Matrix_Convention) is
      A_P : Real_Matrix := (case Convention is
                               when Column_Major => A,
                               when Row_Major => B);
      B_P : Real_Matrix := (case Convention is
                               when Column_Major => B,
                               when Row_Major => A);
      TRANA_P : Real_Trans_Op := (case Convention is
                               when Column_Major => TRANA,
                               when Row_Major => TRANB);
      TRANB_P : Real_Trans_Op := (case Convention is
                               when Column_Major => TRANB,
                               when Row_Major => TRANA);
      M_P : Vector_Size;
      N_P : Vector_Size;
      K_P : Vector_Size;
   begin

      case TRANA_P is
         when No_Transpose =>
            M_P := A_P'Length(2);
            K_P := A_P'Length(1);
         when Transpose =>
            M_P := A_P'Length(1);
            K_P := A_P'Length(2);
      end case;

      case TRANB_P is
         when No_Transpose =>
            N_P := B_P'Length(1);
         when Transpose =>
            N_P := B_P'Length(2);
      end case;

      case Precision is
         when Single =>
            SGEMM(TRANA => Map_Trans_Op(TRANA_P),
                  TRANB => Map_Trans_Op(TRANB_P),
                  M => M_P,
                  N => N_P,
                  K => K_P,
                  ALPHA => ALPHA,
                  A => A_P,
                  LDA => A_P'Length(2),
                  B => B_P,
                  LDB => B_P'Length(2),
                  BETA => BETA,
                  C => C,
                  LDC => C'Length(2)
                  );
         when Double =>
            DGEMM(TRANA => Map_Trans_Op(TRANA_P),
                   TRANB => Map_Trans_Op(TRANB_P),
                  M => M_P,
                  N => N_P,
                  K => K_P,
                  ALPHA => ALPHA,
                  A => A_P,
                  LDA => A_P'Length(2),
                  B => B_P,
                  LDB => B_P'Length(2),
                  BETA => BETA,
                  C => C,
                  LDC => C'Length(2)
                  );
      end case;
   end gemm;

end aBLAS.Real_BLAS;
