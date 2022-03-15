function ang = check_angle(ang)

% Ensures 0 <= ang < 360
% AJ: How about this for code tidy:

ang = mod(ang, 360);

% for i = 1:size(ang1)
%     while( ang1(i) < 0)
%         ang1(i) = ang1(i) + 360;
%     end
%     while ang1(i) >= 360
%         ang1(i) = ang1(i) - 360;
%     end
% end

% ----------------------------------------------------------------------------------------------------------------------
