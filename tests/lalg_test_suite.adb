-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines
-- Unit tests

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

with LALG, LALG.Real_BLAS;

pragma Elaborate_All(LALG, LALG.Real_BLAS);

with BLAS_Real_Level1;
with BLAS_Real_Level2;
with BLAS_Real_Level3;

package body LALG_Test_Suite is
   use AUnit.Test_Suites;

   Result : aliased Test_Suite;

   package Float_BLAS_Base is new LALG(Real => Float);
   package Float_BLAS is new Float_BLAS_Base.Real_BLAS;

   package Float_Real_Level1 is new BLAS_Real_Level1(Float_BLAS_Base,
                                                      Float_BLAS,
                                                      "Float");
   Test_Float_Real_Level1 : aliased Float_Real_Level1.Level1_Test;

   package Float_Real_Level2 is new BLAS_Real_Level2(Float_BLAS_Base,
                                                                 Float_BLAS,
                                                                 "Float");
   Test_Float_Real_Level2 : aliased Float_Real_Level2.Level2_Test;

   package Float_Real_Level3 is new BLAS_Real_Level3(Float_BLAS_Base,
                                                                 Float_BLAS,
                                                                 "Float");
   Test_Float_Real_Level3 : aliased Float_Real_Level3.Level3_Test;

   package Long_Float_BLAS_Base is new LALG(Real => Long_Float);
   package Long_Float_BLAS is new Long_Float_BLAS_Base.Real_BLAS;

   package Long_Float_Real_Level1 is new BLAS_Real_Level1(Long_Float_BLAS_Base,
                                                                      Long_Float_BLAS,
                                                                      "Long_Float");
   Test_Long_Float_Real_Level1 : aliased Long_Float_Real_Level1.Level1_Test;

   package Long_Float_Real_Level2 is new BLAS_Real_Level2(Long_Float_BLAS_Base,
                                                                      Long_Float_BLAS,
                                                                      "Long_Float");
   Test_Long_Float_Real_Level2 : aliased Long_Float_Real_Level2.Level2_Test;

   package Long_Float_Real_Level3 is new BLAS_Real_Level3(Long_Float_BLAS_Base,
                                                                      Long_Float_BLAS,
                                                                      "Long_Float");
   Test_Long_Float_Real_Level3 : aliased Long_Float_Real_Level3.Level3_Test;

   -----------
   -- Suite --
   -----------

   function Suite return AUnit.Test_Suites.Access_Test_Suite is
   begin
      Add_Test (Result'Access, Test_Float_Real_Level1'Access);
      Add_Test (Result'Access, Test_Float_Real_Level2'Access);
      Add_Test (Result'Access, Test_Float_Real_Level3'Access);
      Add_Test (Result'Access, Test_Long_Float_Real_Level1'Access);
      Add_Test (Result'Access, Test_Long_Float_Real_Level2'Access);
      Add_Test (Result'Access, Test_Long_Float_Real_Level3'Access);
      return Result'Access;
   end Suite;

end LALG_Test_Suite;
