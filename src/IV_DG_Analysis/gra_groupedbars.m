function [hB] = gra_groupedbars( M, E, varargin )
% Plots a grouped bar chart, including error bars on bars.
% Grouping of bars follows default behaviour of BAR.
%
%       [barObjHandle] = gra_groupedbars( barHeights, errors );
%       [barObjHandle] = gra_groupedbars( barHeights, errors, axisHandle );


if isempty(varargin) || ~ishandle( varargin{1} )
%     hF = figure;
    ax = axes;
elseif ishandle( varargin{1} )
    ax = varargin{1};
end
hB             = bar(ax, M);
m_shape = size(M);
a = [hB.XOffset];
%X              = bsxfun( @plus, cell2mat(get(hB,'XData').'),  [hB.XOffset] );
offset_array = reshape(repmat(a, m_shape(1), 1), 1, m_shape(1) * m_shape(2));
X              = bsxfun( @plus, cell2mat(get(hB, 'XData').'), offset_array);
X = reshape(X, m_shape(1), m_shape(2));
hold(ax, 'on');
hL             = errorbar( ax, X, M, E, 'k-' );
[hL.LineStyle] = deal( 'none' ); 
end

