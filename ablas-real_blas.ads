-- aBLAS
-- An Ada 2012 binding to BLAS

generic
   type Real is digits <>;
package aBLAS.Real_BLAS is

   type Real_Scalar(Element : access Real'Base) is null record
     with Implicit_Dereference => Element;

   type Real_1D_Array is array (Integer range <>) of aliased Real'Base
     with Pack, Convention => Fortran;
   type Real_2D_Array is array (Integer range <>, Integer range <>) of aliased Real'Base
     with Pack, Convention => Fortran;

   type Real_Vector_Handle is limited private;
   type Real_Vector_Constant_Handle is limited private;
   type Real_Matrix_Handle is limited private;
   type Real_Matrix_Constant_Handle is limited private;

   type Real_Vector is interface;
   function Length(V : Real_Vector) return Positive is abstract;
   function Stride(V : Real_Vector) return Positive is abstract;
   function Handle(V : in out Real_Vector) return Real_Vector_Handle is abstract;
   function Constant_Handle(V : in Real_Vector) return Real_Vector_Constant_Handle is abstract;
   function Item(V : aliased in Real_Vector; I : Integer) return Real is abstract;
   function Variable_Reference(V: aliased in out Real_Vector; I : Integer)
                               return Real_Scalar is abstract;


   type Concrete_Real_Vector(N : Positive) is new Real_Vector with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Length(V : Concrete_Real_Vector) return Positive;
   function Stride(V : Concrete_Real_Vector) return Positive;
   function Handle(V : in out Concrete_Real_Vector) return Real_Vector_Handle;
   function Constant_Handle(V : in Concrete_Real_Vector) return Real_Vector_Constant_Handle;
   function Values(V : Concrete_Real_Vector) return Real_1D_Array;

   function Make(A : Real_1D_Array) return Concrete_Real_Vector;

   type Real_Vector_View(Base : access Concrete_Real_Vector'Class) is new Real_Vector
   with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Length(V : Real_Vector_View) return Positive;
   function Stride(V : Real_Vector_View) return Positive;
   function Handle(V : in out Real_Vector_View) return Real_Vector_Handle;
   function Constant_Handle(V : in Real_Vector_View) return Real_Vector_Constant_Handle;
   function Make(V : access Concrete_Real_Vector'Class;
                 Start : Positive;
                 Stride : Positive;
                 Length : Natural := 0) return Real_Vector_View;

   type Real_Matrix is interface;
   function Rows(V : Real_Matrix) return Positive is abstract;
   function Columns(V : Real_Matrix) return Positive is abstract;
   function Leading_Dimension(V : Real_Matrix) return Positive is abstract;
   function Handle(V : in out Real_Matrix) return Real_Matrix_Handle is abstract;
   function Constant_Handle(V : in Real_Matrix) return Real_Matrix_Constant_Handle is abstract;
   function Item(V : aliased in Real_Matrix; R, C : Integer) return Real is abstract;
   function Variable_Reference(V: aliased in out Real_Matrix; R, C : Integer)
                               return Real_Scalar is abstract;

   type Concrete_Real_Matrix(M, N : Positive) is new Real_Matrix with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Rows(V : Concrete_Real_Matrix) return Positive;
   function Columns(V : Concrete_Real_Matrix) return Positive;
   function Leading_Dimension(V : Concrete_Real_Matrix) return Positive;
   function Handle(V : in out Concrete_Real_Matrix) return Real_Matrix_Handle;
   function Constant_Handle(V : in Concrete_Real_Matrix) return Real_Matrix_Constant_Handle;
   function Item(V : aliased in Concrete_Real_Matrix; R, C : Integer) return Real;
   function Variable_Reference(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                               return Real_Scalar;
   function Make(A : Real_2D_Array) return Concrete_Real_Matrix;

   type Real_Matrix_Vector(Base : access Concrete_Real_Matrix'Class) is new Real_Vector
   with private
     with Constant_Indexing => Item,
     Variable_Indexing => Variable_Reference;
   function Length(V : Real_Matrix_Vector) return Positive;
   function Stride(V : Real_Matrix_Vector) return Positive;
   function Handle(V : in out Real_Matrix_Vector) return Real_Vector_Handle;
   function Constant_Handle(V : in Real_Matrix_Vector) return Real_Vector_Constant_Handle;

   function Row(V : in out Concrete_Real_Matrix'Class; R : Positive) return Real_Matrix_Vector
     with Pre => R <= V.Rows;
   function Column(V : in out Concrete_Real_Matrix'Class; C : Positive) return Real_Matrix_Vector
     with Pre => C <= V.Columns;
   function Trace(V : in out Concrete_Real_Matrix'Class) return Real_Matrix_Vector
     with Pre => V.Rows = V.Columns;

   -- Some equality operators

   function "="(Left : Real_Vector'Class; Right : Real_1D_Array) return Boolean;

   function Approx_Equal(Left : Real_Vector'Class;
                         Right : Real_1D_Array;
                         Epsilon : Real := 0.001) return Boolean;

   function "="(Left : Real_Matrix'Class; Right : Real_2D_Array) return Boolean;

   function Approx_Equal(Left : Real_Matrix'Class;
                         Right : Real_2D_Array;
                         Epsilon : Real := 0.001) return Boolean;

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   -- Generate Givens plane rotation c<-cos(theta), s<-sin(theta) which would
   -- turn a vector [a, b] into [r, 0]. On exit a<-r and b is s or 1/c
   procedure rotg(a, b : in out Real; c, s : out Real);

   subtype Modified_Givens_Params is Real_1D_Array(1..5);

   -- Generate a modified Givens rotation including scaling factors sqrt(d1)
   -- and sqrt(d2). On exit, d1 and d2 are the diagonal elements of the
   -- transformation matrix and x1 is the rotated co-ordinate. 'params'
   -- contains the details necessary to apply the rotation
   procedure rotmg(d1, d2 : in out Real;
                   x1 : in out Real;
                   y1 : in Real;
                   params : out Modified_Givens_Params);

   -- Apply a Givens rotation to X and Y where c=cos(theta) and s=sin(theta)
   procedure rot(X : in out Real_Vector'Class;
                 Y : in out Real_Vector'Class;
                 C : in Real;
                 S : in Real);

   -- Apply a modified Givens rotation to X and Y as specified by the "PARAMS"
   -- generated by rotmg
   procedure rotm(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  PARAMS : in out Modified_Givens_Params);

   -- Y <-> X
   procedure swap(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class);
   -- X <- aX
   procedure scal(X : in out Real_Vector'Class;
                  A : in Real := 1.0);
   -- Y <- X
   procedure copy(X : in Real_Vector'Class;
                  Y : out Real_Vector'Class);

   -- Y <- aX + Y
   procedure axpy(X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  A : in Real := 1.0);

   -- dot <- X^T . Y
   function dot(X, Y : in Real_Vector'Class) return Real;

   -- sdsdot <- X^T . Y + B with accumulation done in extended precision
   -- If Real is already double precision, this is the same as using the regular
   -- dot function and adding B
   function sdsdot(X, Y : in Real_Vector'Class;
                   B : in Real := 0.0) return Real;

   -- nrm2 <- sqrt(X^T . X)
   function nrm2(X : in Real_Vector'Class) return Real;

   --  asum <- |X|_1
   function asum(X : in Real_Vector'Class) return Real;

    --  iamax <- 1st k where X_k = MAX(abs(X_k))
   function iamax(X : in Real_Vector'Class) return Integer;

   -- *************
   -- *************
   -- ** Level 2 **
   -- *************
   -- *************

   -- y <- alpha*A*x + beta*y
   procedure gemv(A : in Real_Matrix'Class;
                  X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANS : in Real_Trans_Op := No_Transpose);

   -- gemv <- alpha*A*x
   function gemv(A : in Real_Matrix'Class;
                 X : in Real_Vector'Class;
                 ALPHA : in Real := 1.0;
                 TRANS : in Real_Trans_Op := No_Transpose)
                 return Real_Vector'Class;

private

   package IntFort renames Interfaces.Fortran;

   Precision : constant Precision_Specification :=
     (if Real'Base'Digits = IntFort.Real'Base'Digits and
        Real'Base'Size = IntFort.Real'Base'Size then Single
      elsif Real'Base'Digits = IntFort.Double_Precision'Base'Digits and
        Real'Base'Size = IntFort.Double_Precision'Base'Size then Double
      else raise Program_Error with "Precision not supported for interfacing with Fortran code");

   -- The Constant_Handle types only give an indication to the Ada compiler of
   -- the mode of the related FORTRAN parameters. There is no way for the Ada
   -- compiler to enforce that the FORTRAN code does not modify data that it
   -- should not.
   type Real_Handle is access all Real'Base;
   type Real_Constant_Handle is access constant Real'Base;

   type Real_Vector_Handle is new Real_Handle;
   type Real_Vector_Constant_Handle is new Real_Constant_Handle;
   type Real_Matrix_Handle is new Real_Handle;
   type Real_Matrix_Constant_Handle is new Real_Constant_Handle;

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
         Stride : Positive;
         Length : Positive;
         Handle : Real_Vector_Handle;
      end record;
   function Item(V : aliased in Real_Vector_View; I : Integer)
                 return Real with Inline;
   function Variable_Reference(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar with Inline;

   type Concrete_Real_Matrix(M, N : Positive) is new Real_Matrix with
      record
         Data : Real_2D_Array(1..M, 1..N);
      end record;

   type Real_Matrix_Vector(Base : access Concrete_Real_Matrix'Class) is new Real_Vector with
      record
         Start_Row, Start_Column : Positive;
         Offset_Row, Offset_Column : Natural;
         Stride : Positive;
         Length : Positive;
         Handle : Real_Vector_Handle;
      end record;

   function Item(V : aliased in Real_Matrix_Vector; I : Integer)
                    return Real with Inline;
   function Variable_Reference(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar with Inline;


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
