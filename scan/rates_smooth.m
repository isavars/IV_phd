function [smMap]=rates_smooth(map,k,varargin)
% Smooth a place map. 
%       [smoothed_map]=rates_smooth(map,k)
%       [smoothed_map]=rates_smooth(map,k,unvis)
% map=map, k=boxcar kernel, unvis=index of unvisited bins (otherwise assume 0=unvis).

if size(map,2)>1
    %%% Place Maps %%%
    if max(size(k))==1;  
        k=ones(k);  % Expand single parameter to flat k-by-k square
    end
    if ~isempty(varargin)
        unvis=varargin{1};
    else
        unvis=map==0;
    end
    map(unvis)=0;
    visTemplate=ones(size(map));
    visTemplate(unvis)=0;
    filtMap=imfilter(map,k);
    filtVis=imfilter(visTemplate,k);
    warning('off', 'MATLAB:divideByZero');
    smMap=filtMap./filtVis;
    warning('on', 'MATLAB:divideByZero');
    smMap(unvis)=nan;
    
elseif size(map,2)==1 
    %%% Direction maps %%%
    f=ones(k,1);
    fMap=imfilter(map,f,'circular');
    smMap=fMap./(sum(f));
    
end