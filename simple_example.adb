

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;
use aBLAS;

procedure Simple_Example is


   package BLAS is new aBLAS.Real_BLAS(Real => Float);

   package BLAS_Util is new BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   package FI_Text_IO is new Ada.Text_IO.Integer_IO(Num => BLAS.Fortran_Integer);
   use FI_Text_IO;

   X : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0, 3.0, 4.0, 5.0));
   X_View : aliased Blas.Real_Vector_View := BLAS.Make(X'Access, 2, 2);

begin

   Put_Line("*** Level 1 ***");

   Put("X => "); Put(X); New_Line;
   Put("X_View => "); Put(X_View); New_Line;

   Put("|X|_1 is via BLAS library is: ");
   Put(BLAS.asum(X)); New_Line;
   Put("|X_View|_1 is via BLAS library is: ");
   Put(BLAS.asum(X_View)); New_Line;

   Put_Line("Scaling X_View by 2.0");
   BLAS.scal(X_View, 2.0);
   Put("X => "); Put(X); New_Line;
   Put("X_View => "); Put(X_View); New_Line;

end Simple_Example;
