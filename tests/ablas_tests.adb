-- aBLAS
-- An Ada 2012 binding to BLAS
-- Unit test runner

with aBLAS_Test_Suite;

with AUnit.Run;
with AUnit.Reporter.Text;

procedure aBLAS_tests is
   procedure Run is new AUnit.Run.Test_Runner (aBLAS_Test_Suite.Suite);
   Reporter : AUnit.Reporter.Text.Text_Reporter;
begin
   Run (Reporter);
end aBLAS_tests;
