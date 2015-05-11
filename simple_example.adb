

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;
use aBLAS;

procedure Simple_Example is


   package BLAS is new aBLAS.Real_BLAS(Real => Float);

   package BLAS_Util is new BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   X : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0, 3.0, 4.0, 5.0));
   X_View1 : aliased Blas.Real_Vector_View := BLAS.Make(X'Access, 2, 2);
   X_View2 : aliased Blas.Real_Vector_View := BLAS.Make(X'Access, 2, 1, 3);

begin

   Put_Line("*** Level 1 ***");

   Put("X => "); Put(X); New_Line;
   Put("X_View1 => "); Put(X_View1); New_Line;
   Put("X_View2 => "); Put(X_View2); New_Line;

   Put("|X|_1 is via BLAS library is: ");
   Put(BLAS.asum(X)); New_Line;
   Put("|X_View1|_1 is via BLAS library is: ");
   Put(BLAS.asum(X_View1)); New_Line;

   Put_Line("Scaling X_View1 by 2.0");
   BLAS.scal(X_View1, 2.0);
   Put("X => "); Put(X); New_Line;
   Put("X_View1 => "); Put(X_View1); New_Line;
   Put("X_View2 => "); Put(X_View2); New_Line;

end Simple_Example;
