-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests


with aBLAS_Test_Suite.Real_Level1;

package body aBLAS_Test_Suite is
   use AUnit.Test_Suites;

   Result : aliased Test_Suite;

   Test_Real_Level1 : aliased ablas_Test_Suite.Real_Level1.Level1_Test;
   -----------
   -- Suite --
   -----------

   function Suite return AUnit.Test_Suites.Access_Test_Suite is
   begin
      Add_Test (Result'Access, Test_Real_Level1'Access);
      return Result'Access;
   end Suite;

end aBLAS_Test_Suite;
