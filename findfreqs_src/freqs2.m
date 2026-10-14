%%%%%%%%%%%%%%%%%%%%%%%%
% user must supply these three values  
%%%%%%%%%%%%%%%%%%%%%%%%%
M_ann = 0.2686; 
Rtip =  12.2221/12; %in ft 
T = 526.38; 
beta = sqrt(1- M_ann^2);
%%%%%%%%%%%%%%%
% want   k_0^2 - beta^2*Xmn^2  not equal to zero
% k_0 is  omega * Rtip/sos = 2 * pi * f * Rtip / c 
%input number of freqencies for each frequency range
%
nfreq0 =8; 
nfreq1 =5; 
nfreq2 =4; 
nfreq3 =3; 
nfreq4 =2; 
nfreq5 =2; 
nfreq6 =2; 
nfreq7 =4; 
nfreq8 =3; 
% find the frequencies then that correspond to the edge frequencies 

load   fort.102;   % output from CompXmn.f code 
subber = fort(:,3).^2*beta^2;     % beta^2* Xmn^2, or k_0^2
gamma = 1.4;
R = 1716; % [ft*lbf/(slug*R)]
fac = sqrt(gamma*R*T)/2/pi/Rtip;  
f = sqrt(subber) * fac;
ftot = f; 


% put the frequencies in order and get rid of mulitples
fsort = unique(f);

% describe the gap between each frequency
fdiff = diff(fsort);

maxfreq = max(fsort);

% because the gaps change so much, various portions of the 
% spectrum have to be handled separately  
%%%%
% this is the part that may need some changes
% 	both the frequency chunks and the stated gap can 
%	be adjusted
%update 11/25/21 this part is now automated using desired number of
%frequencies in each range
%%%%
for i = 200:-1:1
istart = min(find(fsort>1));
iend = max(find(fsort<800));

% get locations in frequency array where gap is bigger than set value
idx = find(fdiff(istart:iend)>i);
keepi0 = istart + idx -1;
if length(keepi0) == nfreq0 | length(keepi0) > nfreq0
    i
    break
end
end
for i = 200:-1:1
istart = min(find(fsort>800));
iend = max(find(fsort<2000));

% get locations in frequency array where gap is bigger than set value
idx = find(fdiff(istart:iend)>i);
keepi = istart + idx -1;
if length(keepi) == nfreq1 | length(keepi) > nfreq1
    i
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>2000));
iend = max(find(fsort<3000));

idx = find(fdiff(istart:iend)>i);
keepi2 = istart + idx -1;
if length(keepi2) == nfreq2 | length(keepi2) > nfreq2
    i
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>3000));
iend = max(find(fsort<4000));

idx = find(fdiff(istart:iend)>i);
keepi3 = istart + idx -1;
if length(keepi3) == nfreq3 | length(keepi3) > nfreq3
    i
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>4000));
iend = max(find(fsort<8000));


idx = find(fdiff(istart:iend)>i);
keepi4 = istart + idx -1;
if length(keepi4) == nfreq4 | length(keepi4) > nfreq4
    i
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>8000));
iend = max(find(fsort<12000));

idx = find(fdiff(istart:iend)>i);
keepi5 = istart + idx -1;
if length(keepi5) == nfreq5 | length(keepi5) > nfreq5
    i
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>12000));
iend = max(find(fsort<16000));

idx = find(fdiff(istart:iend)>i);
keepi6 = istart + idx -1;
if length(keepi6) == nfreq6 | length(keepi6) > nfreq6
    break
end
end
for i = 100:-1:1
istart = min(find(fsort>16000));
iend = max(find(fsort<21000));

idx = find(fdiff(istart:iend)>i);
keepi7 = istart + idx -1;
if length(keepi7) == nfreq7 | length(keepi7) > nfreq7
    break
end
end

for i = 100:-1:1
istart = min(find(fsort>21000));
iend = max(find(fsort<25000));

idx = find(fdiff(istart:iend)>i);
keepi8 = istart + idx -1;
if length(keepi8) == nfreq8 | length(keepi8) > nfreq8
    break
end
end

%put all array locations marking the gaps into one matrix
keeptot = [keepi0; keepi; keepi2; keepi3; keepi4; keepi5; keepi6; keepi7 ; keepi8];


%frequencies we want to use that are away from singularities are 
%	 so use frequencies located in the middle of the gaps
fs = fsort(keeptot+1) -0.5;

% these are the singular frequencies (edge frequencies)
fssing = fsort(keeptot);

% this plot just gives a feel for how the valid frequencies are
%	distributed -- based on how it looks one would change the
%	code
plot(fs,'*')

% write out the frequencies that you want to use and that you 
% don't want to use
save('fs.dat','fs','-ascii');
save('fs.sing','fsort','-ascii');
