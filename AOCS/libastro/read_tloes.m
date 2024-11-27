function y = read_tloes(fname)
%READ_TLOES Reads NORAD-style satellite orbit prediction data files
%   READ_TLOES(FNAME) parses a standard NORAD "two-line orbital
%   element set" (TLOES) data set in file FNAME.  Returns an array
%   of structures, one structure per satellite description contained
%   in the file.
%
%   Each structure contains the following fields:
%    .name            Satellite name
%    .satnum          Satellite number in NORAD SatCat
%    .classification  U (unclassified) or C (classified)
%    .designator      A structure describing International data:
%                         .year, .launch, .piece
%    .epoch           A structure describing the measurement epoch:
%                         .year, .day (including fractional part)
%    .motion          A structure describing change in mean motion:
%                         .dt1 (1st deriv), .dt2 (2nd deriv)
%    .drag            Ballistic coefficient of drag
%    .ephemeris       Ephemeris model, general set to 0
%    .element         Element number
%    .inclination     Inclination, degrees
%    .RAAN            Right ascension of ascending node, degrees
%    .eccentricity    Eccentricity
%    .argp            Argument of perigee
%    .mean_anomaly    In degrees
%    .mean_motion     Revolutions per day
%    .revnum          Satellite revolution number at epoch
%
% For reference, see:
%   http://www.amsat.org
%   http://www.celestrak.com
%
% A web search on "NORAD two-line element sets" produces excellent links as well

% D. Orofino, The MathWorks, MATLAB File Exchange, updated November 2001

fid = open_file(fname);
i=0;
while ~feof(fid),
    [next,msg] = read_tle(fid);
    if ~isempty(msg),
        disp(msg);
        fprintf('Skipping satellite entry.\n\n');
        
    else
        if ~isempty(next),
            i=i+1; y(i)=next;
        end
    end
end


% --------------------------------------------------------
function [fid,msg] = open_file(fname)

fid=fopen(fname,'rt');
if fid==-1,
    error('File not found.');
end

% --------------------------------------------------------
function [y,msg] = read_tle(fid)


y=[];
msg='';

% Read in next 3 non-blank lines
%
for i=1:3,
    [c{i},valid] = read_next_nonempty_line(fid,i);
    if isempty(c{i}), return; end
    c{i}=deblank(c{i});
    if ~valid, 
        msg='Invalid characters found.';
        if i>1,  % can append more info
            msg=[msg ' Satellite name "' c{1} '"'];
        end
        return; 
    end
end


% Parse first line (line "0")
% ----------------

y.name = deblank(c{1});
if length(y.name)>24,
    fprintf('Nonstandard satellite name ("%s") exceeds 24 characters in length.\n\n',y.name);
end


% Parse second line (line "1")
% ----------------

% 1: TLE line number, column 1
% 2: satellite number (NORAD SATCAT catalog number) and classification
%        NNNNNC  
%     number is in columns 3-7
%     classificatin is column 8
%     where C is a classification letter: U=unclassified, S=secret
% 3: International designator
%        YYNNNPPP, YY=last 2 digits of launch year, NNN=launch number of year, PPP=piece of launch
%     columns 10-11, 12-14, and 15-17
% 4: Epoch date
%         YYDDD.DDDDDDD, YY=last 2 digits of epoch year, D=day of year including fractional portion
%      columns 19-20, 21-32
% 5: First time derivative of mean motion, divided by 2, in revolutions per (day^2)
%        columns 34-43
% 6: Second time derivative of mean motion, divided by 6, in revs / (day^3)
%        columns 45-52
% 7: BSTAR drag term (ballistic coefficient), or radiation coefficient
%        columns 54-61
% 8: Ephemeris type (orbital prediction model)
%        column 63
%     1=SGP, 2=SGP4, 3=SDP4, 4=SGP8, 5=SDP8, but is usually set to 0
% 9: Element number and checksum
%           EEEEC  EEEE=element number, C=checksum
%            columns 65-68, 69

line1=c{2};
if length(c{2})~= 69,
    msg = ['Incorrect length of line 1 for satellite "' y.name '".'];
    return
end
% 1: TLE line number
linenum=str2num(line1(1));
if linenum~=1,
    msg = ['Invalid line number encountered on first line of data for satellite "' y.name '".'];
    return
end
% 2: satellite number and classification
%        NNNNNU, U=unclassified 
s=line1(3:7);
y.satnum=str2num(s);
if isempty(y.satnum),
    msg = ['Invalid satellite number encountered for satellite "' y.name '".'];
    return
end
y.classification = line1(8);
% 3: International designator
%        YYNNNPPP, YY=last 2 digits of launch year, NNN=launch number of year, PPP=piece of launch
y.designator.year=str2num(line1(10:11));
y.designator.launch=str2num(line1(12:14));
y.designator.piece=fliplr(deblank(fliplr(deblank(line1(15:17)))));
% 4: Epoch date
%         YYDDD.DDDDDDD, YY=last 2 digits of epoch year, D=day of year including fractional portion
y.epoch.year_trunc = str2num(line1(19:20)); % truncated (2-digit) year
if y.epoch.year_trunc >= 57, century = 1900; else century=2000; end
y.epoch.year=y.epoch.year_trunc + century;
y.epoch.day = str2num(line1(21:32));
% 5: First time derivative of mean motion
y.motion.dt1=str2num(line1(34:43));
% 6: Second time derivative of mean motion
snum=fliplr(deblank(fliplr(deblank(line1(45:50)))));
exp=str2num(line1(51:52));
if isempty(exp), exp=0; end
% fac is a factor that moves the decimal point to start of number,
% and takes into account the additional power-of-ten (sexp)
fac = 10.^(-length(snum) + exp);
y.motion.dt2=str2num(snum)*fac;
% 7: BSTAR drag term / radiation coeff
snum=fliplr(deblank(fliplr(deblank(line1(54:59)))));
exp=str2num(line1(60:61));
if isempty(exp), exp=0; end
% fac is a factor that moves the decimal point to start of number,
% and takes into account the additional power-of-ten (sexp)
%
fac = 10.^(-length(snum) + exp);
y.drag=str2num(snum)*fac;

% 8: Ephemeris type
y.ephemeris=str2num(line1(63));

% 9: Element number and checksum
%           EEEEC  EEEE=element number, C=checksum
y.element=str2num(line1(65:68));
checksum.first=str2num(line1(69));

% Verify first checksum, against two possible check values:
chks=compute_checksum(c{2});
if ~any(checksum.first == chks),
    msg = ['Checksum on first line for satellite "' y.name '" does not match computed value.'];
    return
end


% Parse thirdline (line "2")
% ---------------
% 1: TLE line number, column 1
% 2: satellite number
%        NNNNN, column 3-7
% 3: inclination (degrees)
%        NNN.NNNN, col 9-16
% 4: right ascension of ascending node (RAAN) (degrees)
%        NNN.NNNN, col 18-25
% 5: eccentricity, leading decimal assumed
%        NNNNNNN, col 27-33
% 6: argument of perigee (argp) (degrees)
%        NNN.NNNN, col 35-42
% 7: mean anomaly degrees)
%        NNN.NNNN, col 44-51
% 8: mean motion (revs per day) (NN.NNNNNNN),
%    revolution number at epoch (RRRRR),
%    and checksum (C)
%        NN.NNNNNNNNRRRRRC, col 53-63, 64-68, 69
%     (14 digits after decimal point)

line2=c{3};
if length(c{3})~= 69,
    msg = ['Incorrect length of line 2 for satellite name "' y.name '".'];
    return
end

% 1: TLE line number
linenum=str2num(line2(1));
if linenum ~= 2,
    msg = 'Invalid line number encountered.';
    return
end
% 2: satellite number
%        NNNNN, col 3-7
satnum=str2num(line2(3:7));
if isempty(satnum),
    msg = 'Invalid satellite number encountered.';
    return
end
if satnum ~= y.satnum,
    msg = 'Satellite number changed from line 1 to line 2 of data set.';
    return
end
% 3: inclination (degrees)
%        NNN.NNNN, col 9-16
y.inclination = str2num(line2(9:16));
% 4: right ascension (degrees)
%        NNN.NNNN, col 18-25
y.RAAN = str2num(line2(18:25));

% 5: eccentricity, leading decimal assumed
%        NNNNNNN, col 27-33
snum=fliplr(deblank(fliplr(deblank(line2(27:33)))));
% fac is a factor that moves the decimal point to start of number,
% and takes into account the additional power-of-ten (sexp)
fac = 10.^(-length(snum));
y.eccentricity = str2num(snum)*fac;
% 6: argument of perigee (argp) (degrees)
%        NNN.NNNN, col 35-42
y.argp = str2num(line2(35:42));
% 7: mean anomaly (degrees)
%        NNN.NNNN, col 44-51
y.mean_anomaly = str2num(line2(44:51));
% 8: mean motion (revs per day) (NN.NNNNNNN),
%    revolution number at epoch (RRRRR),
%    and checksum (C)
%        NN.NNNNNNNNRRRRRC, col 53-63, 64-68, 69
%     (14 digits after decimal point)
y.mean_motion=str2num(line2(53:63));
y.revnum = str2num(line2(64:68));
checksum.second=str2num(line2(69));

% Verify second checksum:
chks=compute_checksum(c{3});
if ~any(checksum.second == chks),
    msg = ['Checksum on second line for satellite "' y.name '" does not match computed value.'];
    return
end


% -----------------------------------
function chk = compute_checksum(str)
% Checksum is modulo-10 sum of line content
% Add numbers on line, ignoring all spaces, letters, periods, + signs
% minus signs take the value of "1"
%
% Returns vector of two possible check sums, according to
% two slightly differing algorithms

str=str(1:end-1);  % remove checksum value itself!

nstr=fix(str);
n0=fix('0'); n9=fix('9');
inum=find( (nstr >= n0) & (nstr<=n9));
iminus=find(str=='-');

digits = sum(nstr(inum)-n0);
minuses = length(iminus);  % count '1' for each minus sign
chk1 = mod(digits+minuses,10);

% We compute a second checksum according to a modified rule:
iplus = find(str=='+');
pluses = length(iplus);
chk2 = mod(digits+minuses+2*pluses, 10);

chk = [chk1 chk2];

% -----------------------------------------------
function [next_line,valid] = read_next_nonempty_line(fid,linenum)
% Read next 3 text lines from file

next_line='';
while isempty(next_line) & ~feof(fid),
    next_line=fgetl(fid);
end

if linenum>1,
    valid = has_valid_chars(next_line);
else
    valid = 1; % any chars valid for line 1 (sat name)
end

% -----------------------------------------------
function y = has_valid_chars(str)
% Valid chars in a line include only the following:
%     [0-9 A-Z period space + -]  (brackets not included!)

islett = isletter(str);
isnumber = ((str>='0') & (str<='9'));
ispunct  = (str==' ') | (str=='.') | (str=='+') | (str=='-');
y = all(islett | isnumber | ispunct);

% [EOF] read_tloes.m
