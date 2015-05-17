-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests


with aBLAS_Test_Suite.Real_Level1;

package body aBLAS_Test_Suite is
   use AUnit.Test_Suites;

   Result : aliased Test_Suite;

   package Float_BLAS is new aBLAS.Real_BLAS(Real => Float);
   package Float_Real_Level1 is new aBLAS_Test_Suite.Real_Level1(Float_BLAS,
                                                                 "Float");
   Test_Float_Real_Level1 : aliased Float_Real_Level1.Level1_Test;

   package Long_Float_BLAS is new aBLAS.Real_BLAS(Real => Long_Float);
   package Long_Float_Real_Level1 is new aBLAS_Test_Suite.Real_Level1(Long_Float_BLAS,
                                                                      "Long_Float");
   Test_Long_Float_Real_Level1 : aliased Long_Float_Real_Level1.Level1_Test;
   -----------
   -- Suite --
   -----------

   function Suite return AUnit.Test_Suites.Access_Test_Suite is
   begin
      Add_Test (Result'Access, Test_Float_Real_Level1'Access);
      Add_Test (Result'Access, Test_Long_Float_Real_Level1'Access);
      return Result'Access;
   end Suite;

end aBLAS_Test_Suite;
