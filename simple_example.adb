

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;
use aBLAS;

procedure Simple_Example is


   package BLAS is new aBLAS.Real_BLAS(Real => Float);

   package BLAS_Util is new BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   P : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0, 3.0, 4.0, 5.0));
   P_View1 : aliased Blas.Real_Vector_View := BLAS.Make(V => P'Access,
                                                        Start => 2,
                                                        Stride => 2);
   P_View2 : aliased Blas.Real_Vector_View := BLAS.Make(V => P'Access,
                                                        Start => 2,
                                                        Stride => 1,
                                                        Length => 3);

   A : aliased BLAS.Concrete_Real_Matrix := BLAS.Make(((3.0, 4.0), (5.0, 6.0)));
   X : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0));
   Y : aliased BLAS.Concrete_Real_Vector := BLAS.Make((5.0, 5.0));

begin

   Put_Line("*** Level 1 ***");

   Put("P => "); Put(P); New_Line;
   Put("P_View1 => "); Put(P_View1); New_Line;
   Put("P_View2 => "); Put(P_View2); New_Line;

   Put("|P|_1 is via BLAS library is: ");
   Put(BLAS.asum(P)); New_Line;
   Put("|P_View1|_1 is via BLAS library is: ");
   Put(BLAS.asum(P_View1)); New_Line;

   Put_Line("Scaling P_View1 by 2.0");
   BLAS.scal(P_View1, 2.0);
   Put("P => "); Put(P); New_Line;
   Put("P_View1 => "); Put(P_View1); New_Line;
   Put("P_View2 => "); Put(P_View2); New_Line;

   New_Line;
   Put_Line("*** Level 2 ***");
   Put("A => "); Put(A); New_Line;
   Put("A(2,1) => "); Put(A(2,1)); New_Line;
   Put("X => "); Put(X); New_Line;
   Put("Y => "); Put(Y); New_Line;
   Put_Line("Perform Y <- A*X + Y");
   BLAS.gemv(A, X, Y, ALPHA => 1.0, BETA => 1.0);
   Put("Y => "); Put(Y); New_Line;

end Simple_Example;
