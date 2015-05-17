-- aBLAS
-- An Ada 2012 binding to BLAS

with aBLAS.Real_BLAS.Imports;
with aBLAS.Internal;

package body aBLAS.Real_BLAS is

   package Fortran_Imports is new aBLAS.Real_BLAS.Imports;
   use Fortran_Imports;

   --
   -- Concrete_Real_Vector
   --

   function Length(V : Concrete_Real_Vector) return Positive is
     (V.N);

   function Stride(V : Concrete_Real_Vector) return Positive is (1);

   function Handle(V : in out Concrete_Real_Vector) return Real_Vector_Handle is
     (V.Data(1)'Unchecked_Access);

   function Constant_Handle(V : in Concrete_Real_Vector) return Real_Vector_Constant_Handle is
     (V.Data(1)'Unchecked_Access);

   function Values(V : Concrete_Real_Vector) return Real_1D_Array is (V.Data);

   function Item(V : aliased in Concrete_Real_Vector; I : Integer)
                 return Real is (V.Data(I));

   function Variable_Reference(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar is ((Element => V.Data(I)'Access));

   function Make(A : Real_1D_Array) return Concrete_Real_Vector is
     (Concrete_Real_Vector'(N => A'Length, Data => A));

   --
   -- Real_Vector_View
   --

   function Length(V : Real_Vector_View) return Positive is
     (V.Length);

   function Stride(V : Real_Vector_View) return Positive is
     (V.Stride);

   function Handle(V : in out Real_Vector_View) return Real_Vector_Handle is
     (V.Handle);

   function Constant_Handle(V : in Real_Vector_View) return Real_Vector_Constant_Handle is
     (Real_Vector_Constant_Handle(V.Handle));

   function Item(V : aliased in Real_Vector_View; I : Integer)
                 return Real is
     (V.Base.Data(V.Start + (I-1)*V.Stride));

   function Variable_Reference(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar is
     ((Element => V.Base.Data(V.Start + (I-1)*V.Stride)'Access));

   function Make(V : access Concrete_Real_Vector'Class;
                 Start : Positive;
                 Stride : Positive;
                 Length : Natural := 0) return Real_Vector_View is
      Result : Real_Vector_View(Base => V);
      Max_Length : Positive := 1 + (Positive(V.Length) - Start) / Stride;
   begin
      if Start > Positive(V.Length) then
         raise Constraint_Error
           with "Starting position beyond end of vector: " & Integer'Image(Start);
      elsif Length > Max_Length then
         raise Constraint_Error
           with "Vector view too long for actual vector: " & Integer'Image(Length);
      end if;

      Result.Start := Start;
      Result.Stride := Stride;
      Result.Length := (if Length > 0 then Length else Max_Length);

      Result.Handle := V.Data(Start)'Unchecked_Access;
      return Result;
   end Make;

   --
   -- Concrete_Real_Matrix
   --

   function Rows(V : Concrete_Real_Matrix) return Positive is (V.M);
   function Columns(V : Concrete_Real_Matrix) return Positive is (V.N);
   function Leading_Dimension(V : Concrete_Real_Matrix) return Positive is (V.N);
   function Handle(V : in out Concrete_Real_Matrix) return Real_Matrix_Handle is
     (V.Data(1,1)'Unchecked_Access);
   function Constant_Handle(V : in Concrete_Real_Matrix) return Real_Matrix_Constant_Handle is
     (V.Data(1,1)'Unchecked_Access);
   function Item(V : aliased in Concrete_Real_Matrix; R, C : Integer) return Real is
      (V.Data(R, C));
   function Variable_Reference(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                               return Real_Scalar is
      ((Element => V.Data(R,C)'Access));
   function Make(A : Real_2D_Array) return Concrete_Real_Matrix is
      ((M => A'Length(1), N => A'Length(2), Data => A));

   --
   -- Real_Matrix_Vector
   --
   function Length(V : Real_Matrix_Vector) return Positive is
     (V.Length);

   function Stride(V : Real_Matrix_Vector) return Positive is
     (V.Stride);

   function Handle(V : in out Real_Matrix_Vector) return Real_Vector_Handle is
     (V.Handle);

   function Constant_Handle(V : in Real_Matrix_Vector) return Real_Vector_Constant_Handle is
     (Real_Vector_Constant_Handle(V.Handle));

   function Item(V : aliased in Real_Matrix_Vector; I : Integer)
                 return Real is
     (V.Base.Data(V.Start_Row + (I-1)*V.Offset_Row,
                  V.Start_Column + (I-1)*V.Offset_Column ));

   function Variable_Reference(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar is
     ((Element => V.Base.Data(V.Start_Row + (I-1)*V.Offset_Row,
                  V.Start_Column + (I-1)*V.Offset_Column )'Access));

   function Row(V : in out Concrete_Real_Matrix'Class; R : Positive) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => R,
                          Start_Column => 1,
                          Offset_Row => 0,
                          Offset_Column => 1,
                          Stride => V.M,
                          Length => V.N,
                          Handle => V.Data(R, 1)'Unchecked_Access));

    function Column(V : in out Concrete_Real_Matrix'Class; C : Positive) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => 1,
                          Start_Column => C,
                          Offset_Row => 1,
                          Offset_Column => 0,
                          Stride => 1,
                          Length => V.M,
                          Handle => V.Data(1, C)'Unchecked_Access));

    function Trace(V : in out Concrete_Real_Matrix'Class) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => 1,
                          Start_Column => 1,
                          Offset_Row => 1,
                          Offset_Column => 1,
                          Stride => V.M+1,
                          Length => V.M,
                          Handle => V.Data(1, 1)'Unchecked_Access));

   -- Some equality operators

   function "="(Left : Real_Vector'Class; Right : Real_1D_Array) return Boolean is
     (Left.Length = Right'Length and then
        (for all I in 1..Left.Length => Left.Item(I) = Right(Right'First+I-1)));

   function Approx_Equal(Left : Real_Vector'Class;
                         Right : Real_1D_Array;
                         Epsilon : Real := 0.001) return Boolean is
     (
      Left.Length = Right'Length and then
        (for all I in 1..Left.Length => abs(Left.Item(I)-Right(Right'First+I-1)) <= Epsilon)
     );

   function "="(Left : Real_Matrix'Class; Right : Real_2D_Array) return Boolean is
   begin
      if Left.Rows /= Right'Length(1) or Left.Columns /= Right'Length(2) then
         return False;
      else
         for I in 1..Left.Rows loop
            for J in 1..Left.Columns loop
               if Left.Item(I, J) /= Right(Right'First(1)+I-1, Right'First(2)+J-1) then
                  return False;
               end if;
            end loop;
         end loop;
         return true;
      end if;
   end "=";

   function Approx_Equal(Left : Real_Matrix'Class;
                         Right : Real_2D_Array;
                         Epsilon : Real := 0.001) return Boolean is
   begin
      if Left.Rows /= Right'Length(1) or Left.Columns /= Right'Length(2) then
         return False;
      else
         for I in 1..Left.Rows loop
            for J in 1..Left.Columns loop
               if abs(Left.Item(I, J) - Right(Right'First(1)+I-1, Right'First(2)+J-1)) > Epsilon then
                  return False;
               end if;
            end loop;
         end loop;
         return true;
      end if;
   end Approx_Equal;

   Map_Trans_Op : Internal.Map_Trans_Op_Array renames Internal.Map_Trans_Op;

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





--     ---------
--     -- dot --
--     ---------
--
--     function dot
--       (X, Y : in Real_Vector;
--        INCX : in Increment := 1;
--        INCY : in Increment := 1;
--        N : in Vector_Size := 0)
--        return Real
--     is
--        NN : constant Vector_Size := (if N = 0 then X'Length else N);
--     begin
--        case Precision is
--           when Single => return SDOT(NN,X,INCX,Y,INCY);
--           when Double => return DDOT(NN,X,INCX,Y,INCY);
--        end case;
--     end dot;
--
--     ------------
--     -- sdsdot --
--     ------------
--
--     function sdsdot
--       (SX, SY: in Real_Vector;
--        SB : in Real := 0.0;
--        INCX : in Increment := 1;
--        INCY : in Increment := 1;
--        N : in Vector_Size := 0)
--        return Real
--     is
--        NN : constant Vector_Size := (if N = 0 then SX'Length else N);
--     begin
--        case Precision is
--           when Single => return SDSDOT(NN,SB,SX,INCX,SY,INCY);
--           when Double => return DDOT(NN,SX,INCX,SY,INCY) + SB;
--        end case;
--     end sdsdot;
--
--     ----------
--     -- nrm2 --
--     ----------
--
--     function nrm2
--       (X : in Real_Vector;
--        INCX : in Increment := 1;
--        N : in Vector_Size := 0)
--        return Real
--     is
--        NN : constant Vector_Size := (if N = 0 then X'Length else N);
--     begin
--        case Precision is
--           when Single => return SNRM2(NN,X,INCX);
--           when Double => return DNRM2(NN,X,INCX);
--        end case;
--     end nrm2;
--

--
--     function iamax(X : in Real_Vector;
--                    INCX : in Increment := 1;
--                    N : in Vector_Size := 0) return Integer
--     is
--        NN : constant Vector_Size := (if N = 0 then X'Length else N);
--     begin
--        case Precision is
--           when Single => return Natural(ISAMAX(NN,X,INCX))+X'First-1;
--           when Double => return Natural(IDAMAX(NN,X,INCX))+X'First-1;
--        end case;
--     end iamax;
--

--
--     function gemv(A : in Real_Matrix;
--                   X : in Real_Vector;
--                   ALPHA : in Real := 1.0;
--                   TRANS : in Real_Trans_Op := No_Transpose;
--                   INCX : in Increment := 1;
--                   M, N : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--                   return Real_Vector is
--        Y : Real_Vector(1..A'Length(1));
--        -- As Beta is being set to zero, Y should only be written to and not read
--        -- so it does not matter that it is uninitialised.
--     begin
--        gemv(A => A,
--             X => X,
--             Y => Y,
--             ALPHA => ALPHA,
--             BETA => 0.0,
--             TRANS => TRANS,
--             INCX => INCX,
--             INCY => 1,
--             M => M,
--             N => N,
--             Convention => Convention);
--        return Y;
--     end gemv;
--
--     procedure ger(X : in Real_Vector;
--                   Y : in Real_Vector;
--                   A : in out Real_Matrix;
--                   ALPHA : in Real := 1.0;
--                   INCX : in Increment := 1;
--                   INCY : in Increment := 1;
--                   M, N : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--     is
--        M_P : Vector_Size;
--        N_P : Vector_Size;
--     begin
--        case Convention is
--           when Row_Major =>
--              M_P := (if N = 0 then Y'Length else N);
--              N_P := (if M = 0 then X'Length else M);
--              case Precision is
--                 when Single =>
--                    SGER(M => M_P,
--                         N => N_P,
--                         ALPHA => ALPHA,
--                         X => Y,
--                         INCX => INCY,
--                         Y => X,
--                         INCY => INCX,
--                         A => A,
--                         LDA => A'Length(2));
--                 when Double =>
--                    DGER(M => M_P,
--                         N => N_P,
--                         ALPHA => ALPHA,
--                         X => Y,
--                         INCX => INCY,
--                         Y => X,
--                         INCY => INCX,
--                         A => A,
--                         LDA => A'Length(2));
--              end case;
--           when Column_Major =>
--              M_P := (if M = 0 then X'Length else M);
--              N_P := (if N = 0 then Y'Length else N);
--              case Precision is
--                 when Single =>
--                    SGER(M => M_P,
--                         N => N_P,
--                         ALPHA => ALPHA,
--                         X => X,
--                         INCX => INCX,
--                         Y => Y,
--                         INCY => INCY,
--                         A => A,
--                         LDA => A'Length(2));
--                 when Double =>
--                    DGER(M => M_P,
--                         N => N_P,
--                         ALPHA => ALPHA,
--                         X => X,
--                         INCX => INCX,
--                         Y => Y,
--                         INCY => INCY,
--                         A => A,
--                         LDA => A'Length(2));
--              end case;
--        end case;
--     end ger;
--
--     -- *************
--     -- *************
--     -- ** Level 3 **
--     -- *************
--     -- *************
--
--     procedure gemm(A : in Real_Matrix;
--                    B : in Real_Matrix;
--                    C : in out Real_Matrix;
--                    ALPHA : in Real := 1.0;
--                    BETA : in Real := 0.0;
--                    TRANA : in Real_Trans_Op := No_Transpose;
--                    TRANB : in Real_Trans_Op := No_Transpose;
--                    M, N, K : in Vector_Size := 0;
--                    Convention : in Matrix_Convention := Default_Matrix_Convention) is
--        M_P : Vector_Size;
--        N_P : Vector_Size;
--        K_P : Vector_Size;
--     begin
--        case Convention is
--           when Column_Major =>
--              case TRANA is
--              when No_Transpose =>
--                 M_P := (if M = 0 then A'Length(2) else M);
--                 K_P := (if K = 0 then A'Length(1) else K);
--              when Transpose =>
--                 M_P := (if M = 0 then A'Length(1) else M);
--                 K_P := (if K = 0 then A'Length(2) else K);
--              end case;
--
--              case TRANB is
--              when No_Transpose =>
--                 N_P := (if N = 0 then B'Length(1) else N);
--              when Transpose =>
--                 N_P := (if N = 0 then B'Length(2) else N);
--              end case;
--
--              case Precision is
--              when Single =>
--                 SGEMM(TRANA => Map_Trans_Op(TRANA),
--                       TRANB => Map_Trans_Op(TRANB),
--                       M => M_P,
--                       N => N_P,
--                       K => K_P,
--                       ALPHA => ALPHA,
--                       A => A,
--                       LDA => A'Length(2),
--                       B => B,
--                       LDB => B'Length(2),
--                       BETA => BETA,
--                       C => C,
--                       LDC => C'Length(2)
--                      );
--              when Double =>
--                 DGEMM(TRANA => Map_Trans_Op(TRANA),
--                       TRANB => Map_Trans_Op(TRANB),
--                       M => M_P,
--                       N => N_P,
--                       K => K_P,
--                       ALPHA => ALPHA,
--                       A => A,
--                       LDA => A'Length(2),
--                       B => B,
--                       LDB => B'Length(2),
--                       BETA => BETA,
--                       C => C,
--                       LDC => C'Length(2)
--                      );
--              end case;
--
--           when Row_Major =>
--              case TRANA is
--              when No_Transpose =>
--                 M_P := (if M = 0 then B'Length(2) else M);
--                 K_P := (if K = 0 then B'Length(1) else K);
--              when Transpose =>
--                 M_P := (if M = 0 then B'Length(1) else M);
--                 K_P := (if K = 0 then B'Length(2) else K);
--              end case;
--
--              case TRANB is
--              when No_Transpose =>
--                 N_P := (if N = 0 then A'Length(1) else N);
--              when Transpose =>
--                 N_P := (if N = 0 then A'Length(2) else N);
--              end case;
--
--              case Precision is
--              when Single =>
--                 SGEMM(TRANA => Map_Trans_Op(TRANB),
--                       TRANB => Map_Trans_Op(TRANA),
--                       M => M_P,
--                       N => N_P,
--                       K => K_P,
--                       ALPHA => ALPHA,
--                       A => B,
--                       LDA => B'Length(2),
--                       B => A,
--                       LDB => A'Length(2),
--                       BETA => BETA,
--                       C => C,
--                       LDC => C'Length(2)
--                      );
--              when Double =>
--                 DGEMM(TRANA => Map_Trans_Op(TRANB),
--                       TRANB => Map_Trans_Op(TRANA),
--                       M => M_P,
--                       N => N_P,
--                       K => K_P,
--                       ALPHA => ALPHA,
--                       A => B,
--                       LDA => B'Length(2),
--                       B => A,
--                       LDB => A'Length(2),
--                       BETA => BETA,
--                       C => C,
--                       LDC => C'Length(2)
--                      );
--              end case;
--
--        end case;
--
--     end gemm;
--
--     function gemm(A : in Real_Matrix;
--                   B : in Real_Matrix;
--                   ALPHA : in Real := 1.0;
--                   TRANA : in Real_Trans_Op := No_Transpose;
--                   TRANB : in Real_Trans_Op := No_Transpose;
--                   M, N, K : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--                   return Real_Matrix is
--        C : Real_Matrix(1..A'Length(1), 1..B'Length(2));
--        -- As Beta is being set to zero, C should only be written to and not read
--        -- so it does not matter that it is uninitialised.
--     begin
--        gemm(A => A,
--             B => B,
--             C => C,
--             ALPHA => ALPHA,
--             BETA => 0.0,
--             TRANA => TRANA,
--             TRANB => TRANB,
--             M => M,
--             N => N,
--             K => K,
--             Convention => Convention);
--        return C;
--     end gemm;

end aBLAS.Real_BLAS;
