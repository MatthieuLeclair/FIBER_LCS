function [lda1, lda2, xi1x, xi1y, xi2x, xi2y] = CauchyGreen(m11, m12, m21, m22)
   % The inputs of this matrix are the terms of the flow map gradient tensor,
   % M.  The term m11 is the term in the (1,1) index of this matrix.
    
   % Calculate the terms of the symetric Cauchy-Green deformation tensor which
   % is the product [M]^T * M where we have taken a transpose of the first
   % term.
   c11 = m11.^2 + m21.^2;
   c12 = m11.*m12 + m21.*m22;
   c22 = m12.^2 + m22.^2;

   % Eigenvalues of the Cauchy-Green tensor
   lda1 = (c11+c22)/2 - .5*sqrt( (c11-c22).^2 + 4*c12.^2);
   lda2 = (c11+c22)/2 + .5*sqrt( (c11-c22).^2 + 4*c12.^2);

   % Coordinates of xi1
   xi1x = c12./sqrt(c12.^2 + (lda1-c11).^2);
   xi1y = (lda1-c11)./sqrt(c12.^2 + (lda1-c11).^2);

   % In some cases these terms become degenerate.  This occurs when the
   % denominator is zero.  There are two ways to solve for the terms of the
   % eigenvector.  The first attempt has been to use the terms in the top line
   % of the matrix.  If this has failed we attempt the second line of the
   % matrix.
   f = find( sqrt(c12.^2 + (lda1-c11).^2) < 1.e-5 );
   xi1x(f) = (lda1(f)-c22(f))./sqrt(c12(f).^2 + (lda1(f)-c22(f)).^2);
   xi1y(f) = c12(f)./sqrt(c12(f).^2 + (lda1(f)-c22(f)).^2);

   % Now we account for the potential points where degenerencies occur.  It 
   % is possible for the eigenvectors to not be unique.  In this case they
   % will be degenerate and we set the values of the eigenvectors as follows.
   f = isnan(xi1x) & ~isnan(lda2);
   xi1x(f) = 0;
   xi1y(f) = 1;
   
   % Set xi2 as the direct normal vector to xi1
   xi2x = -xi1y;
   xi2y =  xi1x;
end