-- aBLAS
-- An Ada 2012 binding to BLAS

with Ada.Text_IO;

generic
   Default_Columns : Ada.Text_IO.Count := 5;
   Default_Rows : Ada.Text_IO.Count := 5;
   Default_Fore : Ada.Text_IO.Field := 2;
   Default_Aft : Ada.Text_IO.Field := Real'Digits-1;
   Default_Exp : Ada.Text_IO.Field := 3;
package aBLAS.Real_BLAS.Util is

   use Ada.Text_IO;
   subtype Row_Column_Count is Count range 3..Count'Last;

   procedure Put(File : in File_Type;
                 Item : in Real_Vector'Class;
                 Columns : in Row_Column_Count := Default_Columns;
                 Fore : in Field := Default_Fore;
                 Aft : in Field := Default_Aft;
                 Exp : in Field := Default_Exp);

   procedure Put(Item : in Real_Vector'Class;
                 Columns : in Row_Column_Count := Default_Columns;
                 Fore : in Field := Default_Fore;
                 Aft : in Field := Default_Aft;
                 Exp : in Field := Default_Exp);

   procedure Put(File : in File_Type;
                 Item : in Real_Matrix'Class;
                 Columns : in Row_Column_Count := Default_Columns;
                 Rows : in Row_Column_Count := Default_Rows;
                 Fore : in Field := Default_Fore;
                 Aft : in Field := Default_Aft;
                 Exp : in Field := Default_Exp);

   procedure Put(Item : in Real_Matrix'Class;
                 Columns : in Row_Column_Count := Default_Columns;
                 Rows : in Row_Column_Count := Default_Rows;
                 Fore : in Field := Default_Fore;
                 Aft : in Field := Default_Aft;
                 Exp : in Field := Default_Exp);

end aBLAS.Real_BLAS.Util;
