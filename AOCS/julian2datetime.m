function utc_datetime = julian2datetime(julian_date)
    % Convert Julian Date to UTC datetime
    %
    % Inputs:
    %   julian_date - Julian Date
    %
    % Output:
    %   utc_datetime - [year, month, day, hour, minute, second]

    % Constants
    JD_J2000 = 2451545.0; % Julian Date of J2000.0
    days_per_century = 36525; % Days per Julian century

    % Compute the number of days from J2000.0
    days_from_J2000 = julian_date - JD_J2000;

    % Calculate the integer part of the Julian Date
    J = floor(julian_date + 0.5);

    % Calculate the fraction of the day
    F = (julian_date + 0.5) - J;

    % Calculate the year, month, day
    if J >= 2299161
        % Gregorian calendar
        A = floor((J - 1867216.25) / 36524.25);
        B = J + 1 + A - floor(A / 4);
    else
        % Julian calendar
        B = J;
    end

    C = B + 1524;
    D = floor((C - 122.1) / 365.25);
    E = floor(365.25 * D);
    G = floor((C - E) / 30.6001);

    % Day
    day = C - E + F - floor(30.6001 * G);

    % Month
    if G < 13.5
        month = G - 1;
    else
        month = G - 13;
    end

    % Year
    if month > 2.5
        year = D - 4716;
    else
        year = D - 4715;
    end

    % Calculate the hour, minute, and second
    day_fraction = F * 24;
    hour = floor(day_fraction);
    minute = floor((day_fraction - hour) * 60);
    second = (day_fraction - hour - minute / 60) * 3600;

    % Combine into a single datetime array
    utc_datetime = [year, month, floor(day), hour, minute, second];
end
