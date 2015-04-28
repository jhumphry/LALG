

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;
with Ada.Numerics.Generic_Real_Arrays;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;

procedure Simple_Example is
   package GRA is new Ada.Numerics.Generic_Real_Arrays(Real => Float);

   package BLAS is new aBLAS.Real_BLAS(Real => Float,
                                       Real_Vector => GRA.Real_Vector,
                                       Real_Matrix => GRA.Real_Matrix);

   package BLAS_Util is new BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   SX : GRA.Real_Vector(1..3) := (1.0, 2.0, 3.0);
   SY : GRA.Real_Vector(1..3) := (6.0, 5.0, 4.0);

   X  : GRA.Real_Vector(1..2) := (3.0, 4.0);
   A  : GRA.Real_Matrix(1..2, 1..2) := ((1.0, -1.0),
                                        (0.5, 2.0));

   BV : GRA.Real_Vector(1..6) := (others => 3.5);
   BM : GRA.Real_Matrix(1..6, 1..6) := ((1.0, -1.0, others => 0.7),
                                        (0.5, 2.0, others => 0.8),
                                       others => (0.33, 0.67, others => 0.9));

begin
   Put("SX is: ");
   Put(SX); New_Line;
   Put("SY is: ");
   Put(SY); New_Line;

   Put("Dot Product SX.SY is: ");
   Put(BLAS.Dot(SX, SY));
   New_Line;

   New_Line;
   Put("X is: ");
   Put(X); New_Line;
   Put_Line("A is: ");
   Put(A); New_Line;

   New_Line;
   Put("BV is: ");
   Put(BV);  New_Line;
   Put_Line("BM is: ");
   Put(BM); New_Line;

end Simple_Example;
