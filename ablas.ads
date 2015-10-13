-- aBLAS
-- An Ada 2012 binding to BLAS

with Interfaces.Fortran;

package aBLAS is

   package IntFort renames Interfaces.Fortran;
   use type Interfaces.Fortran.Fortran_Integer;

   -- Fortran precision corresponding to the float type used
   type Precision_Specification is (Single, Double);

   -- Matrix element order convention
   type Matrix_Convention is (Row_Major, Column_Major);
   Ada_Convention : constant Matrix_Convention := Row_Major;
   C_Convention : constant Matrix_Convention := Row_Major;
   Fortran_Convention : constant Matrix_Convention := Column_Major;

   -- Enumeration types for operation specifications

   -- Operand side for non-commutative operations
   type Side is (Left, Right);

   -- Use the upper or lower half of a symmetric matrix
   type UpLo is (Upper, Lower);

   -- Transpose operation specifier
   type Trans_Op is (No_Transpose, Transpose, Conj_Transpose);
   subtype Real_Trans_Op is Trans_Op range No_Transpose..Transpose;

   -- Conjugate operation specifier
   type Conj_Op is (No_Conj, Conj);

   -- Describe diagonal matrix
   type Diag is (Non_Unit_Diag, Unit_Diag);

   -- JRot type
   type JRot is (Inner, Outer, Sorted);

   -- Increments can never be zero
   subtype Incr is IntFort.Fortran_Integer
   with Static_Predicate => Incr in
     IntFort.Fortran_Integer'First..-1|+1..IntFort.Fortran_Integer'Last;

   -- 0 is used as a sentinel for the natural array size
   subtype Vector_Size is IntFort.Fortran_Integer range 0..IntFort.Fortran_Integer'Last;

end aBLAS;
