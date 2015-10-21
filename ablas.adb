-- aBLAS
-- An Ada 2012 binding to BLAS

package body aBLAS is

   --
   -- Concrete_Real_Vector
   --

   function Length(V : Concrete_Real_Vector) return Positive is
     (V.N);

   function Stride(V : Concrete_Real_Vector) return Positive is (1);

   function Handle(V : in out Concrete_Real_Vector) return Real_Vector_Handle is
     (V.Data(1)'Unchecked_Access);

   function Constant_Handle(V : in Concrete_Real_Vector) return Real_Vector_Constant_Handle is
     (V.Data(1)'Unchecked_Access);

   function Values(V : Concrete_Real_Vector) return Real_1D_Array is (V.Data);

   function Item(V : aliased in Concrete_Real_Vector; I : Integer)
                 return Real is (V.Data(I));

   function Variable_Reference(V: aliased in out Concrete_Real_Vector; I : Integer)
                               return Real_Scalar is ((Element => V.Data(I)'Access));

   function Make(A : Real_1D_Array) return Concrete_Real_Vector is
     (Concrete_Real_Vector'(N => A'Length, Data => A));

   function Zeros(Length : Positive) return Concrete_Real_Vector is
     (Concrete_Real_Vector'(N => Length, Data => (others => 0.0)));

   function Ones(Length : Positive) return Concrete_Real_Vector is
     (Concrete_Real_Vector'(N => Length, Data => (others => 1.0)));

   --
   -- Real_Vector_View
   --

   function Length(V : Real_Vector_View) return Positive is
     (V.Length);

   function Stride(V : Real_Vector_View) return Positive is
     (V.Stride);

   function Handle(V : in out Real_Vector_View) return Real_Vector_Handle is
     (V.Handle);

   function Constant_Handle(V : in Real_Vector_View) return Real_Vector_Constant_Handle is
     (Real_Vector_Constant_Handle(V.Handle));

   function Item(V : aliased in Real_Vector_View; I : Integer)
                 return Real is
     (V.Base.Data(V.Start + (I-1)*V.Stride));

   function Variable_Reference(V: aliased in out Real_Vector_View; I : Integer)
                               return Real_Scalar is
     ((Element => V.Base.Data(V.Start + (I-1)*V.Stride)'Access));

   function Make(V : access Concrete_Real_Vector'Class;
                 Start : Positive;
                 Stride : Positive;
                 Length : Natural := 0) return Real_Vector_View is
      Result : Real_Vector_View(Base => V);
      Max_Length : constant Positive := 1 + (V.Length - Start) / Stride;
   begin
      if Start > V.Length then
         raise Constraint_Error
           with "Starting position beyond end of vector: " & Integer'Image(Start);
      elsif Length > Max_Length then
         raise Constraint_Error
           with "Vector view too long for actual vector: " & Integer'Image(Length);
      end if;

      Result.Start := Start;
      Result.Stride := Stride;
      Result.Length := (if Length > 0 then Length else Max_Length);

      Result.Handle := V.Data(Start)'Unchecked_Access;
      return Result;
   end Make;

   --
   -- Concrete_Real_Matrix
   --

   function Rows(V : Concrete_Real_Matrix) return Positive is (V.M);
   function Columns(V : Concrete_Real_Matrix) return Positive is (V.N);
   function Leading_Dimension(V : Concrete_Real_Matrix) return Positive is (V.M);
   function Handle(V : in out Concrete_Real_Matrix) return Real_Matrix_Handle is
     (V.Data(1,1)'Unchecked_Access);
   function Constant_Handle(V : in Concrete_Real_Matrix) return Real_Matrix_Constant_Handle is
     (V.Data(1,1)'Unchecked_Access);
   function Item(V : aliased in Concrete_Real_Matrix; R, C : Integer) return Real is
      (V.Data(R, C));
   function Variable_Reference(V: aliased in out Concrete_Real_Matrix; R, C : Integer)
                               return Real_Scalar is
      ((Element => V.Data(R,C)'Access));
   function Make(A : Real_2D_Array) return Concrete_Real_Matrix is
      ((M => A'Length(1), N => A'Length(2), Data => A));
   function Zeros(Rows, Columns : Positive) return Concrete_Real_Matrix is
     ((M => Rows, N => Columns, Data => (others => (others => 0.0))));
   function Ones(Rows, Columns : Positive) return Concrete_Real_Matrix is
     ((M => Rows, N => Columns, Data => (others => (others => 1.0))));
   function Identity(Rows : Positive) return Concrete_Real_Matrix is
      Result : Concrete_Real_Matrix(Rows, Rows);
   begin
      for I in 1..Rows loop
         for J in 1..Rows loop
            Result.Data(I,J) := (if I=J then 1.0 else 0.0);
         end loop;
      end loop;
      return Result;
   end Identity;

   --
   -- Real_Matrix_Vector
   --
   function Length(V : Real_Matrix_Vector) return Positive is
     (V.Length);

   function Stride(V : Real_Matrix_Vector) return Positive is
     (V.Stride);

   function Handle(V : in out Real_Matrix_Vector) return Real_Vector_Handle is
     (V.Handle);

   function Constant_Handle(V : in Real_Matrix_Vector) return Real_Vector_Constant_Handle is
     (Real_Vector_Constant_Handle(V.Handle));

   function Item(V : aliased in Real_Matrix_Vector; I : Integer)
                 return Real is
     (V.Base.Data(V.Start_Row + (I-1)*V.Offset_Row,
                  V.Start_Column + (I-1)*V.Offset_Column ));

   function Variable_Reference(V: aliased in out Real_Matrix_Vector; I : Integer)
                               return Real_Scalar is
     ((Element => V.Base.Data(V.Start_Row + (I-1)*V.Offset_Row,
                  V.Start_Column + (I-1)*V.Offset_Column )'Access));

   function Row(V : in out Concrete_Real_Matrix'Class; R : Positive) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => R,
                          Start_Column => 1,
                          Offset_Row => 0,
                          Offset_Column => 1,
                          Stride => V.M,
                          Length => V.N,
                          Handle => V.Data(R, 1)'Unchecked_Access));

    function Column(V : in out Concrete_Real_Matrix'Class; C : Positive) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => 1,
                          Start_Column => C,
                          Offset_Row => 1,
                          Offset_Column => 0,
                          Stride => 1,
                          Length => V.M,
                          Handle => V.Data(1, C)'Unchecked_Access));

    function Diagonal(V : in out Concrete_Real_Matrix'Class) return Real_Matrix_Vector is
     (Real_Matrix_Vector'(Base => V'Access,
                          Start_Row => 1,
                          Start_Column => 1,
                          Offset_Row => 1,
                          Offset_Column => 1,
                          Stride => V.M+1,
                          Length => V.M,
                          Handle => V.Data(1, 1)'Unchecked_Access));

   --
   -- Packed_Real_Matrix
   --
   function Rows(V : Packed_Real_Matrix) return Positive is (V.M);
   function Columns(V : Packed_Real_Matrix) return Positive is (V.M);
   function Handle(V : in out Packed_Real_Matrix) return Real_Packed_Matrix_Handle is
     (V.Data(1)'Unchecked_Access);
   function Constant_Handle(V : in Packed_Real_Matrix) return Real_Packed_Matrix_Constant_Handle is
     (V.Data(1)'Unchecked_Access);

   --
   -- Symmetric_Real_Matrix
   --

   function Item(V : aliased in Symmetric_Real_Matrix; R, C : Integer) return Real is
      S : constant Integer := Integer'Min(R, C);
      T : constant Integer := Integer'Max(R, C);
      U : constant Integer := (T * (T-1)) / 2 + S;
   begin
      return V.Data(U);
   end Item;

   function Variable_Reference(V: aliased in out Symmetric_Real_Matrix; R, C : Integer)
                               return Real_Scalar is
      S : constant Integer := Integer'Min(R, C);
      T : constant Integer := Integer'Max(R, C);
      U : constant Integer := (T * (T-1)) / 2 + S;
   begin
      return Real_Scalar'(Element => V.Data(U)'Access);
   end Variable_Reference;

   function Make(A : Real_2D_Array) return Symmetric_Real_Matrix is
      K : Positive := 1;
   begin
      return R : Symmetric_Real_Matrix(M => A'Length(1),
                                       L => (A'Length(1) * (A'Length(1) + 1)) / 2) do
         for I in A'Range(2) loop
            for J in A'First(1)..A'First(1)+(I-A'First(2)) loop
               R.Data(K) := A(J, I);
               K := K + 1;
            end loop;
         end loop;
      end return;
   end Make;

   function Zeros(Rows : Positive) return Symmetric_Real_Matrix is
     (Symmetric_Real_Matrix'(M => Rows,
                             L => (Rows * (Rows + 1))/2,
                             Data => (others => 0.0)));

   function Ones(Rows : Positive) return Symmetric_Real_Matrix is
     (Symmetric_Real_Matrix'(M => Rows,
                             L => (Rows * (Rows + 1))/2,
                             Data => (others => 1.0)));

   function Identity(Rows : Positive) return Symmetric_Real_Matrix is
      Index : Integer := 1;
   begin
      return R : Symmetric_Real_Matrix(M => Rows,
                                       L => (Rows * (Rows + 1))/2) do
         R.Data := (others => 0.0);
         for I in 1..Rows loop
            R.Data(Index) := 1.0;
            Index := Index + I + 1;
         end loop;
      end return;
   end Identity;

   -- Some equality operators

   function "="(Left : Real_Vector'Class; Right : Real_1D_Array) return Boolean is
     (Left.Length = Right'Length and then
        (for all I in 1..Left.Length => Left.Item(I) = Right(Right'First+I-1)));

   function Approx_Equal(Left : Real_Vector'Class;
                         Right : Real_1D_Array;
                         Epsilon : Real := 0.001) return Boolean is
     (
      Left.Length = Right'Length and then
        (for all I in 1..Left.Length => abs(Left.Item(I)-Right(Right'First+I-1)) <= Epsilon)
     );

   function "="(Left : Real_Abstract_Matrix'Class;
                Right : Real_2D_Array) return Boolean is
   begin
      if Left.Rows /= Right'Length(1) or Left.Columns /= Right'Length(2) then
         return False;
      else
         for I in 1..Left.Rows loop
            for J in 1..Left.Columns loop
               if Left.Item(I, J) /= Right(Right'First(1)+I-1, Right'First(2)+J-1) then
                  return False;
               end if;
            end loop;
         end loop;
         return True;
      end if;
   end "=";

   function Approx_Equal(Left : Real_Abstract_Matrix'Class;
                         Right : Real_2D_Array;
                         Epsilon : Real := 0.001) return Boolean is
   begin
      if Left.Rows /= Right'Length(1) or Left.Columns /= Right'Length(2) then
         return False;
      else
         for I in 1..Left.Rows loop
            for J in 1..Left.Columns loop
               if abs(Left.Item(I, J) - Right(Right'First(1)+I-1, Right'First(2)+J-1)) > Epsilon then
                  return False;
               end if;
            end loop;
         end loop;
         return True;
      end if;
   end Approx_Equal;

end aBLAS;
