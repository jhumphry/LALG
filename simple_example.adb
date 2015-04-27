

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;
with Ada.Numerics.Generic_Real_Arrays;

with aBLAS, aBLAS.Real_BLAS;

procedure Simple_Example is
   package GRA is new Ada.Numerics.Generic_Real_Arrays(Real => Float);

   package BLAS is new aBLAS.Real_BLAS(Real => Float,
                                       Real_Vector => GRA.Real_Vector,
                                       Real_Matrix => GRA.Real_Matrix);


   procedure Print_Vector(A : in GRA.Real_Vector; Col : Positive := 3) is
      C : Natural := 0;
   begin
      for I in A'Range loop
            Put(A(I));
         C := C + 1;
         if C mod Col = 0 then
            New_Line;
         end if;
      end loop;
      if C mod Col /= 0 then
         New_Line;
      end if;
   end Print_Vector;

   SX : GRA.Real_Vector(1..3) := (1.0, 2.0, 3.0);
   SY : GRA.Real_Vector(1..3) := (4.0, 5.0, 6.0);
   Dot_Product : Float;

begin
   Put_Line("SX is:");
   Print_Vector(SX);
   New_Line;

   Put_Line("SY is:");
   Print_Vector(SY);
   New_Line;

   Dot_Product := BLAS.Dot(SX, SY);

   Put("Dot Product is: ");
   Put(Dot_Product);
   New_Line;

end Simple_Example;
