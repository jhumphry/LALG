-- aBLAS
-- An Ada 2012 binding to BLAS

with Interfaces.Fortran;

private package aBLAS.Internal is

   package IntFort renames Interfaces.Fortran;

   function TF(Item : in Character) return IntFort.Character_Set
               renames IntFort.To_Fortran;

   type Map_Trans_Op_Array is array (Trans_Op) of IntFort.Character_Set;
   Map_Trans_Op : constant Map_Trans_Op_Array := (No_Trans => TF('N'),
                                                  Trans => TF('T'),
                                                  Conj_Trans => TF('C'));

end aBLAS.Internal;
