% Coordinates
BS = [0, 0, 15];
RIS = [35, 20, 15];
BD1 = [30, 5, 0];
BD2 = [40, 5, 0];

% Function to compute azimuth and elevation
ang = @(v) [atan2(v(2), v(1)), atan2(v(3), sqrt(v(1)^2 + v(2)^2))];

% 1. AoD from BS to RIS
v1 = RIS - BS;
[a1, e1] = deal(ang(v1));
fprintf('AoD BS->RIS: Az = %.2f°, El = %.2f°\n', rad2deg(a1), rad2deg(e1));

% 2. AoA at RIS from BS (same vector)
fprintf('AoA RIS<-BS: Az = %.2f°, El = %.2f°\n', rad2deg(a1), rad2deg(e1));

% 3. AoD from RIS to BD1
v2 = BD1 - RIS;
[a2, e2] = deal(ang(v2));
fprintf('AoD RIS->BD1: Az = %.2f°, El = %.2f°\n', rad2deg(a2), rad2deg(e2));

% 4. AoD from RIS to BD2
v3 = BD2 - RIS;
[a3, e3] = deal(ang(v3));
fprintf('AoD RIS->BD2: Az = %.2f°, El = %.2f°\n', rad2deg(a3), rad2deg(e3));