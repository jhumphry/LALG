-- aBLAS
-- An Ada 2012 binding to BLAS

generic
package aBLAS.Real_BLAS is

   -- Enumeration types for operation specifications

   -- Operand side for non-commutative operations
   type Side_Op is (Left, Right);

   -- Use the upper or lower part of a symmetric matrix
   type UpLo_Part is (Upper, Lower);

   -- Transpose operation specifier
   type Trans_Op is (No_Transpose, Transpose, Conj_Transpose);
   subtype Real_Trans_Op is Trans_Op range No_Transpose..Transpose;

   -- Conjugate operation specifier
   type Conj_Op is (No_Conj, Conj);

   -- Describe whether diagonal can be assumed to be all unity
   type Diag_Unit is (Non_Unit_Diag, Unit_Diag);

   -- JRot type
   type JRot is (Inner, Outer, Sorted);

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   -- Generate Givens plane rotation c<-cos(theta), s<-sin(theta) which would
   -- turn a vector [a, b] into [r, 0]. On exit a<-r and b is s or 1/c
   procedure rotg(a, b : in out Real; c, s : out Real);

   subtype Modified_Givens_Params is Real_1D_Array(1..5);

   -- Generate a modified Givens rotation including scaling factors sqrt(d1)
   -- and sqrt(d2). On exit, d1 and d2 are the diagonal elements of the
   -- transformation matrix and x1 is the rotated co-ordinate. 'params'
   -- contains the details necessary to apply the rotation
   procedure rotmg(d1, d2 : in out Real;
                   x1 : in out Real;
                   y1 : in Real;
                   params : out Modified_Givens_Params);

   -- Apply a Givens rotation to X and Y where c=cos(theta) and s=sin(theta)
   procedure rot(X : in out Real_Vector'Class;
                 Y : in out Real_Vector'Class;
                 C : in Real;
                 S : in Real);

   -- Apply a modified Givens rotation to X and Y as specified by the "PARAMS"
   -- generated by rotmg
   procedure rotm(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  PARAMS : in out Modified_Givens_Params);

   -- Y <-> X
   procedure swap(X : in out Real_Vector'Class;
                  Y : in out Real_Vector'Class);
   -- X <- aX
   procedure scal(X : in out Real_Vector'Class;
                  A : in Real := 1.0);
   -- Y <- X
   procedure copy(X : in Real_Vector'Class;
                  Y : out Real_Vector'Class);

   -- Y <- aX + Y
   procedure axpy(X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  A : in Real := 1.0);

   -- dot <- X^T . Y
   function dot(X, Y : in Real_Vector'Class) return Real;

   -- sdsdot <- X^T . Y + B with accumulation done in extended precision
   -- If Real is already double precision, this is the same as using the regular
   -- dot function and adding B
   function sdsdot(X, Y : in Real_Vector'Class;
                   B : in Real := 0.0) return Real;

   -- nrm2 <- sqrt(X^T . X)
   function nrm2(X : in Real_Vector'Class) return Real;

   --  asum <- |X|_1
   function asum(X : in Real_Vector'Class) return Real;

    --  iamax <- 1st k where X_k = MAX(abs(X_k))
   function iamax(X : in Real_Vector'Class) return Integer;

   -- *************
   -- *************
   -- ** Level 2 **
   -- *************
   -- *************

   -- y <- alpha*A*x + beta*y
   procedure gemv(A : in Real_Matrix'Class;
                  X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANS : in Real_Trans_Op := No_Transpose)
     with Pre => (X.Length = A.Columns and Y.Length = A.Rows);

   -- gemv <- alpha*A*x
   function gemv(A : in Real_Matrix'Class;
                 X : in Real_Vector'Class;
                 ALPHA : in Real := 1.0;
                 TRANS : in Real_Trans_Op := No_Transpose)
                 return Real_Vector'Class
     with Pre => (X.Length = A.Columns);

   -- y <- alpha*A*x + beta*y using only the upper or lower triangular part of A
   procedure symv(A : in Real_Matrix'Class;
                  UPLO : in UpLo_Part;
                  X : in Real_Vector'Class;
                  Y : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0)
     with Pre => (A.Columns = A. Rows and
                    X.Length = A.Columns and
                      Y.Length = A.Rows);

   -- symv <- alpha*A*x using only the upper or lower triangular part of A
   function symv(A : in Real_Matrix'Class;
                 UPLO : in UpLo_Part;
                 X : in Real_Vector'Class;
                 ALPHA : in Real := 1.0)
                 return Real_Vector'Class
     with Pre => (A.Columns = A. Rows and
                    X.Length = A.Columns);

   -- A <- alpha*x*yT + A
   procedure ger(X : in Real_Vector'Class;
                 Y : in Real_Vector'Class;
                 A : in out Real_Matrix'Class;
                 ALPHA : in Real := 1.0)
     with Pre => (X.Length = A.Rows and Y.Length = A.Columns);

   -- A <- alpha*x*xT + A using only the upper or lower triangular part of A
   procedure syr(X : in Real_Vector'Class;
                 A : in out Real_Matrix'Class;
                 UPLO : in UpLo_Part;
                 ALPHA : in Real := 1.0)
     with Pre => (A.Rows = A.Columns and
                    X.Length = A.Rows);

   -- A <- alpha*x*yT + alpha*y*xT + A
   -- using only the upper or lower triangular part of A
   procedure syr2(X : in Real_Vector'Class;
                  Y : in Real_Vector'Class;
                  A : in out Real_Matrix'Class;
                  UPLO : in UpLo_Part;
                  ALPHA : in Real := 1.0)
     with Pre => (A.Rows = A.Columns and
                    X.Length = A.Rows and
                      Y.Length = A.Rows);

   -- *************
   -- *************
   -- ** Level 3 **
   -- *************
   -- *************

   -- C <- alpha*A*B + beta*C
   procedure gemm(A : in Real_Matrix'Class;
                  B : in Real_Matrix'Class;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANA : in Real_Trans_Op := No_Transpose;
                  TRANB : in Real_Trans_Op := No_Transpose)
     with Inline;

   -- gemm <- alpha*A*B
   function gemm(A : in Real_Matrix'Class;
                 B : in Real_Matrix'Class;
                 ALPHA : in Real := 1.0;
                 TRANA : in Real_Trans_Op := No_Transpose;
                 TRANB : in Real_Trans_Op := No_Transpose)
                 return Concrete_Real_Matrix
     with Inline;

   -- C <- alpha*A*B + beta*C or C <- alpha*B*A + beta*C, A = AT
   procedure symm(A : in Real_Matrix'Class;
                  SIDE : in Side_Op;
                  UPLO : in UpLo_Part;
                  B : in Real_Matrix'Class;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0)
     with Inline;

   -- symm <- alpha*A*B or symm <- alpha*B*A, A = AT
   function symm(A : in Real_Matrix'Class;
                 SIDE : in Side_Op;
                 UPLO : in UpLo_Part;
                 B : in Real_Matrix'Class;
                 ALPHA : in Real := 1.0)
                 return Concrete_Real_Matrix
     with Inline;

   -- C <- alpha*A*AT + beta*C or C <- alpha*AT*A + beta*C, C=CT
   procedure syrk(A : in Real_Matrix'Class;
                  TRANS : in Real_Trans_Op;
                  UPLO : in UpLo_Part;
                  C : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0)
     with Inline, Pre => (C.Rows = C.Columns and
                            (case TRANS is
                                when No_Transpose => A.Rows = C.Rows,
                                when Transpose => A.Columns = C.Rows)
                         );

   -- syrk <- alpha*A*AT or syrk <- alpha*AT*A
   function syrk(A : in Real_Matrix'Class;
                 TRANS : in Real_Trans_Op;
                 UPLO : in UpLo_Part;
                 ALPHA : in Real := 1.0)
                 return Concrete_Real_Matrix
     with Inline;

   -- C <- alpha*A*BT + alpha*B*AT + beta*C or
   -- C <- alpha*AT*B + alpha*BT*A + beta*C, C=CT
   procedure syr2k(A : in Real_Matrix'Class;
                   B : in Real_Matrix'Class;
                   TRANS : in Real_Trans_Op;
                   UPLO : in UpLo_Part;
                   C : in out Real_Matrix'Class;
                   ALPHA : in Real := 1.0;
                   BETA : in Real := 0.0)
     with Inline, Pre => (C.Rows = C.Columns and
                            (case TRANS is
                                when No_Transpose => A.Rows = B.Rows and
                                  A.Rows = C.Rows,
                                when Transpose => A.Columns = B.Columns and
                                  A.Columns = C.Rows)
                         );

   -- syr2k <- alpha*A*BT + alpha*B*AT or
   -- syr2k <- alpha*AT*B + alpha*BT*A
   function syr2k(A : in Real_Matrix'Class;
                  B : in Real_Matrix'Class;
                  TRANS : in Real_Trans_Op;
                  UPLO : in UpLo_Part;
                  ALPHA : in Real := 1.0)
                  return Concrete_Real_Matrix
     with Inline, Pre => (A.Rows = B.Rows and A.Columns = B.Columns);

   -- B <- alpha*TRANSA(A)*B or B <- alpha*B*TRANSA(A)
   -- A is an upper or lower triangular matrix with unit or non-unit diagonal
   procedure trmm(A : in Real_Matrix'Class;
                  SIDE : in Side_Op;
                  UPLO : in UpLo_Part;
                  TRANSA : in Trans_Op;
                  DIAG : in Diag_Unit;
                  B : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0)
     with Inline;

   -- B <- solve(TRANSA(A) * X = alpha * B) or
   -- B <- solve(X * TRANSA(A) = alpha * B) or
   -- A is an upper or lower triangular matrix with unit or non-unit diagonal
   -- B is a matrix
   procedure trsm(A : in Real_Matrix'Class;
                  SIDE : in Side_Op;
                  UPLO : in UpLo_Part;
                  TRANSA : in Trans_Op;
                  DIAG : in Diag_Unit;
                  B : in out Real_Matrix'Class;
                  ALPHA : in Real := 1.0)
     with Inline;

   -- B <- solve(TRANSA(A) * X = alpha * B) or
   -- B <- solve(X * TRANSA(A) = alpha * B) or
   -- A is an upper or lower triangular matrix with unit or non-unit diagonal
   -- B is a vector
   procedure trsm(A : in Real_Matrix'Class;
                  SIDE : in Side_Op;
                  UPLO : in UpLo_Part;
                  TRANSA : in Trans_Op;
                  DIAG : in Diag_Unit;
                  B : in out Real_Vector'Class;
                  ALPHA : in Real := 1.0)
     with Inline;

end aBLAS.Real_BLAS;
