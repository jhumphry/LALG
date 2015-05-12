-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit tests

with aBLAS, aBLAS.Real_BLAS;
with AUnit.Test_Suites;

package aBLAS_Test_Suite is

   function Suite return AUnit.Test_Suites.Access_Test_Suite;

   package Float_BLAS is new aBLAS.Real_BLAS(Real => Float);

end aBLAS_Test_Suite;
