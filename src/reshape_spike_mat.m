    
    %function needed later from tara's code 
    function [final_mat] = reshape_spike_mat(interp_mat)
    
    interp_mat = permute(interp_mat,[1 3 2]);
    final_mat = reshape(interp_mat,[],51); % Reshape to form 2D matrix
    
    end  