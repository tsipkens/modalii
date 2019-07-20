
clear;
close all;
clc;

% indexFile = '..\Data\Fe (In House, Round 4, Aug)\index.xlsx';
% saveFile = '+inHouse\Fe_signal_vec.mat';

% indexFile = '..\Data\Ag (In House)\index.xlsx';
% saveFile = '+inHouse\Ag_signal_vec_wrise.mat';

indexFile = '..\Data\Mo (In House)\index.xlsx';
saveFile = '+inHouse\Mo_signal_vec.mat';

[~,~,c] = xlsread(indexFile);
c = c(3:end,:);
pathstr = fileparts(indexFile);
addpath(pathstr);

c_sz = size(c);
ind = 1;%1:c_sz(1); % chose whether to selectively load signals
% c = c(ind,:);

signal_vec = [];
for ii=ind
    fname = c{ii,6};
    l = eval(c{ii,3});
    signal = inHouse.loadFiles_ih(fname,l);
    signal.F0 = c{ii,4};
    signal.matl = c{ii,1};
    signal.gas = c{ii,2};
    signal.type = 'Artium\inHouse';
    signal_vec = [signal_vec;signal];
end

disp(' ');
disp('Saving data...');
save(saveFile,'signal_vec');
disp('Save complete.');

