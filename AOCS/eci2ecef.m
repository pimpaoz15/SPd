function r_ecef = eci2ecef(r_eci, julian_date)
    % ECI to ECEF Conversion
    %
    % Inputs:
    %   r_eci       - (3x1) Vector in ECI frame
    %   julian_date - Julian date
    %
    % Output:
    %   r_ecef      - (3x1) Vector in ECEF frame

    % Constants
    JD_2000 = 2451545.0; % Julian Date of J2000.0
    sec_per_day = 86400; % Seconds per day

    % Compute the number of days from J2000.0
    days_from_J2000 = julian_date - JD_2000;

    % Compute the Greenwich Mean Sidereal Time (GMST) in seconds
    GMST_sec = 67310.54841 + (876600 * 3600 + 8640184.812866) * (days_from_J2000 / 36525) + 0.093104 * (days_from_J2000 / 36525)^2 - 6.2e-6 * (days_from_J2000 / 36525)^3;

    % Normalize GMST to be within 0 to 86400 seconds (1 day)
    GMST_sec = mod(GMST_sec, sec_per_day);

    % Convert GMST from seconds to radians
    GMST_rad = (GMST_sec / sec_per_day) * 2 * pi;

    % Rotation matrix from ECI to ECEF
    R = [cos(GMST_rad), sin(GMST_rad), 0;
         -sin(GMST_rad), cos(GMST_rad), 0;
         0, 0, 1];

    % Convert ECI coordinates to ECEF coordinates
    r_ecef = R * r_eci(:);
end
