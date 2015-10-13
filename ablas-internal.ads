-- aBLAS
-- An Ada 2012 binding to BLAS

with Interfaces.Fortran;

private generic package aBLAS.Internal is

   package IntFort renames Interfaces.Fortran;

   function TF(Item : in Character) return IntFort.Character_Set
               renames IntFort.To_Fortran;

   type Map_Trans_Op_Array is array (Trans_Op) of IntFort.Character_Set;
   Map_Trans_Op : constant Map_Trans_Op_Array := (No_Transpose => TF('N'),
                                                  Transpose => TF('T'),
                                                  Conj_Transpose => TF('C'));

end aBLAS.Internal;
