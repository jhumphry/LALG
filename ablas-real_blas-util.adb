-- aBLAS
-- An Ada 2012 binding to BLAS

with Ada.Strings.Fixed;

package body aBLAS.Real_BLAS.Util is

   package Real_Text_IO is new Ada.Text_IO.Float_IO(Num => Real);
   use Real_Text_IO;

   ---------
   -- Put --
   ---------

   procedure Put
     (File : in File_Type;
      Item : in Real_Vector'Class;
      Columns : in Row_Column_Count := Default_Columns;
      Fore : in Field := Default_Fore;
      Aft : in Field := Default_Aft;
      Exp : in Field := Default_Exp)
   is
      N : constant Integer := Item.Length;
   begin
      Put(File, "[");
      if N <= Integer(Columns) then
         for I in 1..N loop
            Put(File => File,
                Item => Item.Item(I),
                Fore => Fore,
                Aft => Aft,
                Exp => Exp);
            Put(File, " ");
         end loop;
      else
         for I in 1..Integer(Columns)-2 loop
            Put(File => File,
                Item => Item.Item(I),
                Fore => Fore,
                Aft => Aft,
                Exp => Exp);
            Put(File, " ");
         end loop;
         Put(File, Ada.Strings.Fixed."*"(Fore+Aft+Exp+1,"."));
         Put(File, " ");
         Put(File => File,
             Item => Item.Item(N),
             Fore => Fore,
             Aft => Aft,
             Exp => Exp);
      end if;
      Put(File, "]");
   end Put;

   ---------
   -- Put --
   ---------

   procedure Put
     (Item : in Real_Vector'Class;
      Columns : in Row_Column_Count := Default_Columns;
      Fore : in Field := Default_Fore;
      Aft : in Field := Default_Aft;
      Exp : in Field := Default_Exp)
   is
   begin
      Put(Standard_Output, Item, Columns, Fore, Aft, Exp);
   end Put;

   ---------
   -- Put --
   ---------

   procedure Put_Matrix_Row
     (File : in File_Type;
      Item : in Abstract_Real_Matrix'Class;
      Row : in Integer;
      Columns : in Row_Column_Count;
      Fore : in Field;
      Aft : in Field;
      Exp : in Field)
   is
      N : constant Integer := Item.Columns;
   begin
      Put(File, "[");
      if N <= Integer(Columns) then
         for I in 1..Item.Columns loop
            Put(File => File,
                Item => Item.Item(Row, I),
                Fore => Fore,
                Aft => Aft,
                Exp => Exp);
            Put(File, " ");
         end loop;
      else
         for I in 1..Integer(Columns)-2 loop
            Put(File => File,
                Item => Item.Item(Row, I),
                Fore => Fore,
                Aft => Aft,
                Exp => Exp);
            Put(File, " ");
         end loop;
         Put(File, Ada.Strings.Fixed."*"(Fore+Aft+Exp+1,"."));
         Put(File, " ");
         Put(File => File,
             Item => Item.Item(Row, Item.Columns),
             Fore => Fore,
             Aft => Aft,
             Exp => Exp);
      end if;
      Put(File, "]");
   end Put_Matrix_Row;

   procedure Put
     (File : in File_Type;
      Item : in Abstract_Real_Matrix'Class;
      Columns : in Row_Column_Count := Default_Columns;
      Rows : in Row_Column_Count := Default_Rows;
      Fore : in Field := Default_Fore;
      Aft : in Field := Default_Aft;
      Exp : in Field := Default_Exp)
   is
      N : constant Integer := Item.Rows;
   begin
      Set_Col(File, 1);
      Put(File, "[");
      if N <= Integer(Rows) then
         for I in 1..Item.Rows loop
            Set_Col(File, 2);
               Put_Matrix_Row(File => File,
                   Item => Item,
                   Row => I,
                   Columns => Columns,
                   Fore => Fore,
                   Aft => Aft,
                   Exp => Exp);
         end loop;
      else
         for I in 1..Integer(Rows)-2 loop
            Set_Col(File, 2);
               Put_Matrix_Row(File => File,
                   Item => Item,
                   Row => I,
                   Columns => Columns,
                   Fore => Fore,
                   Aft => Aft,
                   Exp => Exp);
         end loop;
         Set_Col(File, 3);
         Put(File, Ada.Strings.Fixed."*"(Fore+Aft+Exp+1,"."));
         New_Line(File);
         Set_Col(File, 2);
         Put_Matrix_Row(File => File,
             Item => Item,
             Row => Item.Rows,
             Columns => Columns,
             Fore => Fore,
             Aft => Aft,
             Exp => Exp);
      end if;
      Put(File, "]");
   end Put;

   ---------
   -- Put --
   ---------

   procedure Put
     (Item : in Abstract_Real_Matrix'Class;
      Columns : in Row_Column_Count := Default_Columns;
      Rows : in Row_Column_Count := Default_Rows;
      Fore : in Field := Default_Fore;
      Aft : in Field := Default_Aft;
      Exp : in Field := Default_Exp)
   is
   begin
      Put(Standard_Output, Item, Columns, Rows, Fore, Aft, Exp);
   end Put;

end aBLAS.Real_BLAS.Util;
