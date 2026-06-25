function designMatrix = getDesignMatrixFromBasisFunctions(designMatrixBase,basisFunctions,nShift)
    designMatrix = zeros(size(designMatrixBase,1)+size(basisFunctions,1)-1,size(designMatrixBase,2)*size(basisFunctions,2));
    for i = 1:size(basisFunctions,2)
        idx = i+(0:size(designMatrixBase,2)-1)*size(basisFunctions,2);
        designMatrix(:,idx) = conv2(designMatrixBase,basisFunctions(:,i),'full');
    end
    designMatrix = circshift(designMatrix,-nShift);
end