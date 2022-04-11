-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines
-- Unit test runner

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

with LALG_Test_Suite;

with AUnit.Run;
with AUnit.Reporter.Text;

procedure LALG_tests is
   procedure Run is new AUnit.Run.Test_Runner (LALG_Test_Suite.Suite);
   Reporter : AUnit.Reporter.Text.Text_Reporter;
begin
   Run (Reporter);
end LALG_tests;
