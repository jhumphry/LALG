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
     (SX, SY: in Real_Vector;
      INCX : in Increment := 1;
      INCY : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then SX'Length else N);
   begin
      case Precision is
         when Single => return SDOT(NN,SX,INCX,SY,INCY);
         when Double => return DDOT(NN,SX,INCX,SY,INCY);
         when Unsupported => raise Program_Error;
      end case;
   end dot;

   ----------
   -- nrm2 --
   ----------

   function nrm2
     (SX : in Real_Vector;
      INCX : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : constant Vector_Size := (if N = 0 then SX'Length else N);
   begin
      case Precision is
         when Single => return SNRM2(NN,SX,INCX);
         when Double => return DNRM2(NN,SX,INCX);
         when Unsupported => raise Program_Error;
      end case;
   end nrm2;

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
         when Unsupported => raise Program_Error;
      end case;
   end gemv;


end aBLAS.Real_BLAS;
