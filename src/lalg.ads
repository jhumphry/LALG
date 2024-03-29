-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

private with Interfaces.Fortran;

generic
   type Real is digits <>;
package LALG is

   -- Use the upper or lower part of a matrix
   -- Store the upper of lower part of a packed symmetrical matrix
   type UpLo_Part is (Upper, Lower);

   -- Basic scalar, vector and matrix types

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
   type Packed_Real_Matrix_Handle is limited private;
   type Packed_Real_Matrix_Constant_Handle is limited private;

   type Real_Vector is interface;
   function Length(V : Real_Vector) return Positive is abstract;
   function Stride(V : Real_Vector) return Positive is abstract;
   function Handle(V : in out Real_Vector) return Real_Vector_Handle is abstract;
   function Constant_Handle(V : in Real_Vector) return Real_Vector_Constant_Handle is abstract;
   function Item(V : aliased in Real_Vector; I : Integer) return Real is abstract
     with Pre'Class => (I <= V.Length);
   function Variable_Reference(V: aliased in out Real_Vector; I : Integer)
                               return Real_Scalar is abstract
     with Pre'Class => (I <= V.Length);

   type Concrete_Real_Vector(N : Positive) is new Real_Vector with private
     with Constant_Indexing => Item_CRV,
     Variable_Indexing => Variable_Reference_CRV;
   function Length(V : Concrete_Real_Vector) return Positive;
   function Stride(V : Concrete_Real_Vector) return Positive;
   function Handle(V : in out Concrete_Real_Vector) return Real_Vector_Handle;
   function Constant_Handle(V : in Concrete_Real_Vector) return Real_Vector_Constant_Handle;
   function Values(V : Concrete_Real_Vector) return Real_1D_Array;
   function Make(A : Real_1D_Array) return Concrete_Real_Vector;
   function Zeros(Length : Positive) return Concrete_Real_Vector;
   function Ones(Length : Positive) return Concrete_Real_Vector;
   function Item_CRV(V : aliased in Concrete_Real_Vector; I : Integer) return Real;
   function Variable_Reference_CRV(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar;

   type Real_Vector_View(Base : access Concrete_Real_Vector'Class) is new Real_Vector
   with private
     with Constant_Indexing => Item_RVV,
     Variable_Indexing => Variable_Reference_RVV;
   function Length(V : Real_Vector_View) return Positive;
   function Stride(V : Real_Vector_View) return Positive;
   function Handle(V : in out Real_Vector_View) return Real_Vector_Handle;
   function Constant_Handle(V : in Real_Vector_View) return Real_Vector_Constant_Handle;
   function Make(V : access Concrete_Real_Vector'Class;
                 Start : Positive;
                 Stride : Positive;
                 Length : Natural := 0) return Real_Vector_View;
   function Item_RVV(V : aliased in Real_Vector_View; I : Integer)
                 return Real;
   function Variable_Reference_RVV(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar;

   type Abstract_Real_Matrix is interface;
   function Rows(V : Abstract_Real_Matrix) return Positive is abstract;
   function Columns(V : Abstract_Real_Matrix) return Positive is abstract;
   function Item(V : aliased in Abstract_Real_Matrix; R, C : Integer) return Real is abstract
     with Pre'Class => (R <= V.Rows and C <= V.Columns);
   function Variable_Reference(V: aliased in out Abstract_Real_Matrix; R, C : Integer)
                               return Real_Scalar is abstract
     with Pre'Class => (R <= V.Rows and C <= V.Columns);

   type Real_Matrix is interface and Abstract_Real_Matrix;
   function Rows(V : Real_Matrix) return Positive is abstract;
   function Columns(V : Real_Matrix) return Positive is abstract;
   function Leading_Dimension(V : Real_Matrix) return Positive is abstract;
   function Handle(V : in out Real_Matrix) return Real_Matrix_Handle is abstract;
   function Constant_Handle(V : in Real_Matrix) return Real_Matrix_Constant_Handle is abstract;

   type Concrete_Real_Matrix(M, N : Positive) is new Real_Matrix with private
     with Constant_Indexing => Item_CRM,
     Variable_Indexing => Variable_Reference_CRM;
   function Rows(V : Concrete_Real_Matrix) return Positive;
   function Columns(V : Concrete_Real_Matrix) return Positive;
   function Leading_Dimension(V : Concrete_Real_Matrix) return Positive;
   function Handle(V : in out Concrete_Real_Matrix) return Real_Matrix_Handle;
   function Constant_Handle(V : in Concrete_Real_Matrix) return Real_Matrix_Constant_Handle;
   function Make(A : Real_2D_Array) return Concrete_Real_Matrix;
   function Zeros(Rows, Columns : Positive) return Concrete_Real_Matrix;
   function Ones(Rows, Columns : Positive) return Concrete_Real_Matrix;
   function Identity(Rows : Positive) return Concrete_Real_Matrix;
   function Item_CRM(V : aliased in Concrete_Real_Matrix; R, C : Integer)
                     return Real;
   function Variable_Reference_CRM(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                                return Real_Scalar;

   type Real_Matrix_Vector(Base : access Concrete_Real_Matrix'Class) is new Real_Vector
   with private
     with Constant_Indexing => Item_RMV,
     Variable_Indexing => Variable_Reference_RMV;
   function Length(V : Real_Matrix_Vector) return Positive;
   function Stride(V : Real_Matrix_Vector) return Positive;
   function Handle(V : in out Real_Matrix_Vector) return Real_Vector_Handle;
   function Constant_Handle(V : in Real_Matrix_Vector) return Real_Vector_Constant_Handle;

   function Row(V : in out Concrete_Real_Matrix'Class; R : Positive) return Real_Matrix_Vector
     with Pre => R <= V.Rows;
   function Column(V : in out Concrete_Real_Matrix'Class; C : Positive) return Real_Matrix_Vector
     with Pre => C <= V.Columns;
   function Diagonal(V : in out Concrete_Real_Matrix'Class) return Real_Matrix_Vector
     with Pre => V.Rows = V.Columns;
   function Item_RMV(V : aliased in Real_Matrix_Vector; I : Integer)
                 return Real;
   function Variable_Reference_RMV(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar;

   type Packed_Real_Matrix(<>) is abstract tagged private;
   function Rows(V : Packed_Real_Matrix) return Positive;
   function Columns(V : Packed_Real_Matrix) return Positive;
   function Handle(V : in out Packed_Real_Matrix) return Packed_Real_Matrix_Handle;
   function Constant_Handle(V : in Packed_Real_Matrix) return Packed_Real_Matrix_Constant_Handle;

   type Symmetric_Real_Matrix is new Packed_Real_Matrix and Abstract_Real_Matrix with private
     with Constant_Indexing => Item_SRM,
     Variable_Indexing => Variable_Reference_SRM;
   function Item(V : aliased in Symmetric_Real_Matrix; R, C : Integer)
                 return Real with Inline;
   function Variable_Reference(V : aliased in out Symmetric_Real_Matrix; R, C : Integer)
                               return Real_Scalar with Inline;
   function Make(A : Real_2D_Array; UpLo : UpLo_Part) return Symmetric_Real_Matrix
     with Pre => (A'Length(1) = A'Length(2));
   function Zeros(Rows : Positive; UpLo : UpLo_Part) return Symmetric_Real_Matrix;
   function Ones(Rows : Positive; UpLo : UpLo_Part) return Symmetric_Real_Matrix;
   function Identity(Rows : Positive; UpLo : UpLo_Part) return Symmetric_Real_Matrix;
   function Item_SRM(V : aliased in Symmetric_Real_Matrix; R, C : Integer)
                     return Real;
   function Variable_Reference_SRM(V: aliased in out Symmetric_Real_Matrix; R, C : Integer)
                                return Real_Scalar;

   type Triangular_Real_Matrix is new Packed_Real_Matrix and Abstract_Real_Matrix with private
     with Constant_Indexing => Item_TRM,
     Variable_Indexing => Variable_Reference_TRM;
   function Item(V : aliased in Triangular_Real_Matrix; R, C : Integer)
                 return Real with Inline;
   function Variable_Reference(V : aliased in out Triangular_Real_Matrix; R, C : Integer)
                               return Real_Scalar with Inline;
   function Make(A : Real_2D_Array; UpLo : UpLo_Part) return Triangular_Real_Matrix
     with Pre => (A'Length(1) = A'Length(2));
   function Zeros(Rows : Positive; UpLo : UpLo_Part) return Triangular_Real_Matrix;
   function Ones(Rows : Positive; UpLo : UpLo_Part) return Triangular_Real_Matrix;
   function Identity(Rows : Positive; UpLo : UpLo_Part) return Triangular_Real_Matrix;
   function Item_TRM(V : aliased in Triangular_Real_Matrix; R, C : Integer)
                     return Real;
   function Variable_Reference_TRM(V: aliased in out Triangular_Real_Matrix; R, C : Integer)
                                return Real_Scalar;

   -- Some equality operators

   function "="(Left : Real_Vector'Class; Right : Real_1D_Array) return Boolean;

   function Approx_Equal(Left : Real_Vector'Class;
                         Right : Real_1D_Array;
                         Epsilon : Real := 0.001) return Boolean;

   function "="(Left : Abstract_Real_Matrix'Class;
                Right : Real_2D_Array) return Boolean;

   function Approx_Equal(Left : Abstract_Real_Matrix'Class;
                         Right : Real_2D_Array;
                         Epsilon : Real := 0.001) return Boolean;

private

   package IntFort renames Interfaces.Fortran;

   -- Fortran precision corresponding to the float type used
   type Precision_Specification is (Single, Double);

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
   type Packed_Real_Matrix_Handle is new Real_Handle;
   type Packed_Real_Matrix_Constant_Handle is new Real_Constant_Handle;

   type Concrete_Real_Vector(N : Positive) is new Real_Vector with
      record
         Data : Real_1D_Array(1..N);
      end record;

   function Item(V : aliased in Concrete_Real_Vector; I : Integer)
                 return Real with Inline;
   function Item_CRV(V : aliased in Concrete_Real_Vector; I : Integer)
                     return Real renames Item;
   function Variable_Reference(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar with Inline;
   function Variable_Reference_CRV(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar renames Variable_Reference;

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
   function Item_RVV(V : aliased in Real_Vector_View; I : Integer)
                 return Real renames Item;
   function Variable_Reference(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar with Inline;
   function Variable_Reference_RVV(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar renames Variable_Reference;

   type Concrete_Real_Matrix(M, N : Positive) is new Real_Matrix with
      record
         Data : Real_2D_Array(1..M, 1..N);
      end record;

   function Item(V : aliased in Concrete_Real_Matrix; R, C : Integer) return Real;
   function Item_CRM(V : aliased in Concrete_Real_Matrix; R, C : Integer)
                     return Real renames Item;
   function Variable_Reference(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                               return Real_Scalar;
   function Variable_Reference_CRM(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                                return Real_Scalar renames Variable_Reference;

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
   function Item_RMV(V : aliased in Real_Matrix_Vector; I : Integer)
                 return Real renames Item;
   function Variable_Reference(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar with Inline;
   function Variable_Reference_RMV(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar renames Variable_Reference;

   type Packed_Real_Matrix(M, L : Positive) is tagged
      record
         UpLo : UpLo_Part;
         Data : Real_1D_Array(1..L);
      end record;

   type Symmetric_Real_Matrix is new Packed_Real_Matrix and Abstract_Real_Matrix
   with null record;

   function Item_SRM(V : aliased in Symmetric_Real_Matrix; R, C : Integer)
                     return Real renames Item;
   function Variable_Reference_SRM(V: aliased in out Symmetric_Real_Matrix; R, C : Integer)
                                return Real_Scalar renames Variable_Reference;

   type Triangular_Real_Matrix is new Packed_Real_Matrix and Abstract_Real_Matrix
   with null record;

   function Item_TRM(V : aliased in Triangular_Real_Matrix; R, C : Integer)
                     return Real renames Item;
   function Variable_Reference_TRM(V: aliased in out Triangular_Real_Matrix; R, C : Integer)
                                return Real_Scalar renames Variable_Reference;

end LALG;
