-- aBLAS
-- An Ada 2012 binding to BLAS

with Interfaces.Fortran;
private with Interfaces.C.Pointers;

generic
   type Real is digits <>;
package aBLAS.Real_BLAS is

   subtype Fortran_Integer is Interfaces.Fortran.Fortran_Integer;

   type Real_Scalar(Element : access Real'Base) is null record
     with Implicit_Dereference => Element;

   type Real_1D_Array is array (Integer range <>) of aliased Real'Base
     with Pack, Convention => Fortran;

   type Real_Vector_Handle is limited private;

   type Real_Vector is interface;
   function Length(V : Real_Vector) return Fortran_Integer is abstract;
   function Stride(V : Real_Vector) return Fortran_Integer is abstract;
   function Handle(V : in out Real_Vector) return Real_Vector_Handle is abstract;
   function Item(V : aliased in Real_Vector; I : Integer) return Real is abstract;
   function Variable_Reference(V: aliased in out Real_Vector; I : Integer)
                               return Real_Scalar is abstract;

   type Concrete_Real_Vector(N : Positive) is new Real_Vector with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Length(V : Concrete_Real_Vector) return Fortran_Integer;
   function Stride(V : Concrete_Real_Vector) return Fortran_Integer;
   function Handle(V : in out Concrete_Real_Vector) return Real_Vector_Handle;
   function Values(V : Concrete_Real_Vector) return Real_1D_Array;

   function Make(A : Real_1D_Array) return Concrete_Real_Vector;

   type Real_Vector_View(Base : access Concrete_Real_Vector'Class) is new Real_Vector
   with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Length(V : Real_Vector_View) return Fortran_Integer;
   function Stride(V : Real_Vector_View) return Fortran_Integer;
   function Handle(V : in out Real_Vector_View) return Real_Vector_Handle;

   function Make(V : access Concrete_Real_Vector'Class;
                 Start : Positive;
                 Stride : Positive;
                 Length : Natural := 0) return Real_Vector_View;

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   subtype Modified_Givens_Params is Real_1D_Array(1..5);

   --  asum <- |X|_1
   function asum(X : in out Real_Vector'Class) return Real;

   -- X <- aX
   procedure scal(X : in out Real_Vector'Class;
                  A : in Real := 1.0);

private

   package IntFort renames Interfaces.Fortran;

   Precision : constant Precision_Specification :=
     (if Real'Base'Digits = IntFort.Real'Base'Digits and
        Real'Base'Size = IntFort.Real'Base'Size then Single
      elsif Real'Base'Digits = IntFort.Double_Precision'Base'Digits and
        Real'Base'Size = IntFort.Double_Precision'Base'Size then Double
      else raise Program_Error with "Precision not supported for interfacing with Fortran code");

   package RA_Ptrs is new Interfaces.C.Pointers(Index => Integer,
                                                Element => Real'Base,
                                                Element_Array => Real_1D_Array,
                                                Default_Terminator => 0.0);
   type Real_Vector_Handle is new RA_Ptrs.Pointer;

   type Concrete_Real_Vector(N : Positive) is new Real_Vector with
      record
         Data : Real_1D_Array(1..N);
      end record;
   function Item(V : aliased in Concrete_Real_Vector; I : Integer)
                               return Real with Inline;
   function Variable_Reference(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar with Inline;


   type Real_Vector_View(Base : access Concrete_Real_Vector'Class) is new Real_Vector
   with
      record
         Start : Positive;
         Stride : Integer;
         Length : Positive;
         Handle : Real_Vector_Handle;
      end record;
   function Item(V : aliased in Real_Vector_View; I : Integer)
                 return Real with Inline;
   function Variable_Reference(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar with Inline;

   --     -- *************
   --     -- *************
   --     -- ** Level 1 **
   --     -- *************
--     -- *************
--
--     -- Generate Givens plane rotation c<-cos(theta), s<-sin(theta) which would
--     -- turn a vector [a, b] into [r, 0]. On exit a<-r and b is s or 1/c
--     procedure rotg(a, b : in out Real; c, s : out Real)
--       with Inline;
--
--     subtype Modified_Givens_Params is Real_Vector(1..5);
--
--     -- Generate a modified Givens rotation including scaling factors sqrt(d1)
--     -- and sqrt(d2). On exit, d1 and d2 are the diagonal elements of the
--     -- transformation matrix and x1 is the rotated co-ordinate. 'params'
--     -- contains the details necessary to apply the rotation
--     procedure rotmg(d1, d2 : in out Real;
--                     x1 : in out Real;
--                     y1 : in Real;
--                     params : out Modified_Givens_Params)
--       with Inline;
--
--     -- Apply a Givens rotation to X and Y where c=cos(theta) and s=sin(theta)
--     procedure rot(X : in out Real_Vector;
--                   Y : in out Real_Vector;
--                   C : in Real;
--                   S : in Real;
--                   INCX : in Increment := 1;
--                   INCY : in Increment := 1;
--                   N : in Vector_Size := 0)
--       with Inline;
--
--     -- Apply a modified Givens rotation to X and Y as specified by the "PARAMS"
--     -- generated by rotmg
--     procedure rotm(X : in out Real_Vector;
--                    Y : in out Real_Vector;
--                    PARAMS : in out Modified_Givens_Params;
--                    INCX : in Increment := 1;
--                    INCY : in Increment := 1;
--                    N : in Vector_Size := 0)
--       with Inline;
--
--     -- Y <-> X
--     procedure swap(X : in out Real_Vector;
--                    Y : in out Real_Vector;
--                    INCX : in Increment := 1;
--                    INCY : in Increment := 1;
--                    N : in Vector_Size := 0)
--       with Inline;
--

--
--     -- Y <- X
--     procedure copy(X : in Real_Vector;
--                    Y : out Real_Vector;
--                    INCX : in Increment := 1;
--                    INCY : in Increment := 1;
--                    N : in Vector_Size := 0)
--       with Inline;
--
--     -- Y <- aX + Y
--     procedure axpy(X : in Real_Vector;
--                    Y : in out Real_Vector;
--                    A : in Real := 1.0;
--                    INCX : in Increment := 1;
--                    INCY : in Increment := 1;
--                    N : in Vector_Size := 0)
--       with Inline;
--
--     -- dot <- X^T . Y
--     function dot(X, Y : in Real_Vector;
--                  INCX : in Increment := 1;
--                  INCY : in Increment := 1;
--                  N : in Vector_Size := 0) return Real
--       with Inline;
--
--     -- sdsdot <- SX^T . SY + SB with accumulation done in extended precision
--     -- If Real is already double precision, this is the same as using the regular
--     -- dot function and adding SB
--     function sdsdot(SX, SY : in Real_Vector;
--                     SB : in Real := 0.0;
--                     INCX : in Increment := 1;
--                     INCY : in Increment := 1;
--                     N : in Vector_Size := 0) return Real
--       with Inline;
--
--     -- nrm2 <- sqrt(X^T . X)
--     function nrm2(X : in Real_Vector;
--                   INCX : in Increment := 1;
--                   N : in Vector_Size := 0) return Real
--       with Inline;
--
--     --  asum <- |X|_1
--     function asum(X : in Real_Vector;
--                   INCX : in Increment := 1;
--                   N : in Vector_Size := 0) return Real
--       with Inline;
--
--     --  iamax <- 1st k where X_k = MAX(abs(X_k))
--     function iamax(X : in Real_Vector;
--                    INCX : in Increment := 1;
--                    N : in Vector_Size := 0) return Integer
--       with Inline;
--
--     -- *************
--     -- *************
--     -- ** Level 2 **
--     -- *************
--     -- *************
--
--     -- y <- alpha*A*x + beta*y
--     procedure gemv(A : in Real_Matrix;
--                    X : in Real_Vector;
--                    Y : in out Real_Vector;
--                    ALPHA : in Real := 1.0;
--                    BETA : in Real := 0.0;
--                    TRANS : in Real_Trans_Op := No_Transpose;
--                    INCX, INCY : in Increment := 1;
--                    M, N : in Vector_Size := 0;
--                    Convention : in Matrix_Convention := Default_Matrix_Convention)
--       with Inline;
--
--     -- gemv <- alpha*A*x
--     function gemv(A : in Real_Matrix;
--                   X : in Real_Vector;
--                   ALPHA : in Real := 1.0;
--                   TRANS : in Real_Trans_Op := No_Transpose;
--                   INCX : in Increment := 1;
--                   M, N : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--                   return Real_Vector
--       with Inline;
--
--     -- A <- alpha*x*yT + A
--     procedure ger(X : in Real_Vector;
--                   Y : in Real_Vector;
--                   A : in out Real_Matrix;
--                   ALPHA : in Real := 1.0;
--                   INCX : in Increment := 1;
--                   INCY : in Increment := 1;
--                   M, N : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--       with Inline;
--
--     -- *************
--     -- *************
--     -- ** Level 3 **
--     -- *************
--     -- *************
--
--     -- C <- alpha*A*B + beta*C
--     procedure gemm(A : in Real_Matrix;
--                    B : in Real_Matrix;
--                    C : in out Real_Matrix;
--                    ALPHA : in Real := 1.0;
--                    BETA : in Real := 0.0;
--                    TRANA : in Real_Trans_Op := No_Transpose;
--                    TRANB : in Real_Trans_Op := No_Transpose;
--                    M, N, K : in Vector_Size := 0;
--                    Convention : in Matrix_Convention := Default_Matrix_Convention)
--       with Inline;
--
--     -- gemm <- alpha*A*B
--     function gemm(A : in Real_Matrix;
--                   B : in Real_Matrix;
--                   ALPHA : in Real := 1.0;
--                   TRANA : in Real_Trans_Op := No_Transpose;
--                   TRANB : in Real_Trans_Op := No_Transpose;
--                   M, N, K : in Vector_Size := 0;
--                   Convention : in Matrix_Convention := Default_Matrix_Convention)
--                   return Real_Matrix
--       with Inline;

end aBLAS.Real_BLAS;
