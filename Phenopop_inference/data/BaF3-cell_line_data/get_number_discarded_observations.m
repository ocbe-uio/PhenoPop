function SUM=get_number_discarded_observations(data)

% discarded observations is a function which counts the number of
% observations used by vectorized_inference_k_subpopulations_two_noise_levels. 
% Entire replicates are discarded if the first time point is NaN. Otherwise
% individual observations are discarded if they are NaN. 
[NR NC NT]=size(data);

%get number of discarded observations
stencil=sum(isnan(data(:,:,1)), 'all');
SUM=0;
for t = 1:NT
    canvas = sum(isnan(data(:,:,t)), 'all');    %number of NaN at time t
    if  canvas == stencil
        SUM=SUM+stencil;
        
    elseif canvas > stencil
        newnan= canvas-stencil;
        SUM=SUM + stencil + newnan;
    end
end