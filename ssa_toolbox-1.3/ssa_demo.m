function ssa_demo
%SSA_DEMO        SSA demonstration on synthetic data.
%
%Shows an example of SSA on synthetic data with one stationary and one non-stationary source.
%The upper panel shows the two input signals. Note that due to the mixing, both input signals
%appear non-stationary. The time series is divided into four epochs, for which you can see 
%scatter plots in the middle panel. Here, the coordinate system is chosen according to the 
%decomposition found by SSA. The vertical axis corresponds to the estimated stationary source
%and the horizontal axis to the estimated non-stationary source. You can see that both the 
%mean and the variance changes on the horizontal axis whereas it remains constant on the vertical
%axis. The bottom panel shows the estimated stationary and non-stationary sources.
%
%author 
%  saputra@cs.tu-berlin.de
%
%license
%  This software is distributed under the BSD license. See COPYING for
%  details.

% Copyright (c) 2010, Jan Saputra Müller, Paul von Bünau, Frank C. Meinecke,
% Franz J. Kiraly and Klaus-Robert Müller.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice, this
% list of conditions and the following disclaimer in the documentation and/or other
%  materials provided with the distribution.
% 
% * Neither the name of the Berlin Institute of Technology (Technische Universität
% Berlin) nor the names of its contributors may be used to endorse or promote
% products derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
%  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% fixed random seed
rand('seed', 23);
randn('seed', 34);

%%%%%%%%%%%%%%%%%%%%
% generate data set
%%%%%%%%%%%%%%%%%%%%
n = 500; % number of samples
n_changepoints = 5; % number of changepoints

% random change points
changepoints = sort(floor(n * rand(n_changepoints, 1)));

% generate stationary and non-stationary signal
X = zeros(2, n);
for i=1:(n_changepoints + 1)
    if i == 1
        X(:, 1:changepoints(1)) = mvnrnd([0, 3*randn], randcov, changepoints(1))';
    elseif i == (n_changepoints + 1)
        X(:, (changepoints(end) + 1):end) = mvnrnd([0, 3*randn], randcov, n - changepoints(end))';
    else
        X(:, (changepoints(i-1) + 1):changepoints(i)) = mvnrnd([0, 3*randn], randcov, changepoints(i) - changepoints(i-1))';
    end
end

% random mixing matrix
% random angles in [(1/8)*pi, (3/8)*pi]
alpha = (pi/4)*rand + (1/8)*pi;
beta = (pi/4)*rand + (1/8)*pi;
A = [cos(alpha) sin(alpha); cos(beta) sin(beta)] .* sign(randn(2,2));

% mix signals
Xmixed = A * X;

%%%%%%%%%%%%%%%%%%%%
% run SSA
%%%%%%%%%%%%%%%%%%%%
% 1 stationary source, 5 repetitions, 4 equally-sized epochs, fixed random seed
[est_Ps, est_Pn, est_As, est_An, ssa_results] = ssa(Xmixed, 1, 'reps', 5, 'equal_epochs', 4, 'random_seed', 12345);

% project to stationary and non-stationary subspace
Xest_s = est_Ps * Xmixed; 
Xest_ns = est_Pn * Xmixed;

%%%%%%%%%%%%%%%%%%%%
% plot signals
%%%%%%%%%%%%%%%%%%%%
figure;
sp_ho = 0.1;
sp_ve = 0.05;
s_he = 0.18;
s_sc = 0.02;
s_sc_bo = 0.03;
text_h = 0.025;
text_pos = 0.001;

markers = (n/4):(n/4):((3/4)*n);
% first signal
axes('position', [sp_ho, 1 - sp_ve - s_he, 1 - 2*sp_ho, s_he]);
plot(Xmixed(1, :));
set(gca, 'XTick', [], 'YTick', []);
ylabel('Signal 1');
line_markers(markers);

% second signal
axes('position', [sp_ho, 1 - sp_ve - 2*s_he, 1 - 2*sp_ho, s_he]);
plot(Xmixed(2, :));
set(gca, 'XTick', [], 'YTick', []);
ylabel('Signal 2');
line_markers(markers);

% Plot epoch labels.
for i=1:4  
  loc = [ sp_ho + (i-1)*(1-2*sp_ho)/4, 1-sp_ve+text_pos, (1-2*sp_ho)/4, text_h ];
  annotation(gcf, 'textbox','String',{ sprintf('Epoch %d', i) }, 'Position', loc, 'LineStyle', 'none', 'HorizontalAlignment', 'center');
end

% scatter plots of the epochs
vw = [-3.5 3.5];
eps = n / 4; % epoch size
cbas = 'black'; % color for basis of s- and n-space
cebas = 'red'; % color for estimated basis of s- and n-space
yl = 'Estimated stat. src.';
xl = 'Estimated non-stat. src.';

axes('position', [sp_ho + s_sc, 1 - sp_ve - 3*s_he + s_sc_bo, (1 - 2*sp_ho)/4 - 2*s_sc, s_he - 2*s_sc]);
scatter(Xest_ns(1:eps), Xest_s(1:eps), '.');
box on;
axis equal;
axis square;
set(gca, 'Xlim', vw, 'Ylim', vw, 'XTick', [], 'YTick', []);
xlabel(xl); ylabel(yl);

axes('position', [sp_ho + s_sc + (1 - 2*sp_ho)/4, 1 - sp_ve - 3*s_he + s_sc_bo, (1 - 2*sp_ho)/4 - 2*s_sc, s_he - 2*s_sc]);
scatter(Xest_ns((eps + 1):(2*eps)), Xest_s((eps + 1):(2*eps)), '.');
box on;
axis equal;
axis square;
set(gca, 'Xlim', vw, 'Ylim', vw, 'XTick', [], 'YTick', []);
xlabel(xl); ylabel(yl);

axes('position', [sp_ho + s_sc + 2*(1 - 2*sp_ho)/4, 1 - sp_ve - 3*s_he + s_sc_bo, (1 - 2*sp_ho)/4 - 2*s_sc, s_he - 2*s_sc]);
scatter(Xest_ns((2*eps + 1):(3*eps)), Xest_s((2*eps + 1):(3*eps)), '.');
box on;
axis equal;
axis square;
set(gca, 'Xlim', vw, 'Ylim', vw, 'XTick', [], 'YTick', []);
xlabel(xl); ylabel(yl);

axes('position', [sp_ho + s_sc + 3*(1 - 2*sp_ho)/4, 1 - sp_ve - 3*s_he + s_sc_bo, (1 - 2*sp_ho)/4 - 2*s_sc, s_he - 2*s_sc]);
scatter(Xest_ns((3*eps + 1):end), Xest_s((3*eps + 1):end), '.');
box on;
axis equal;
axis square;
set(gca, 'Xlim', vw, 'Ylim', vw, 'XTick', [], 'YTick', []);
xlabel(xl); ylabel(yl);

set(gcf, 'Position', [20 20 880 980]);

% estimated stationary and non-stationary signal
axes('position', [sp_ho, 1 - sp_ve - 4*s_he, 1 - 2*sp_ho, s_he]);
plot(Xest_s);
set(gca, 'XTick', [], 'YTick', []);
ylabel('Estimated stationary source');

axes('position', [sp_ho, 1 - sp_ve - 5*s_he, 1 - 2*sp_ho, s_he]);
plot(Xest_ns);
set(gca, 'XTick', [], 'YTick', []);
ylabel('Estimated non-stationary source');

% random covariance matrix with a fixed stationary part
function C = randcov
sig = 1;
sc = 1;
m = sqrt(sig*sig/2);
A = [m m; sc*randn(1, 2)];
C = A*A';

function drawaxis(x, y, c, s)
h1 = line([0, s*x(1)], [0, s*x(2)]);
h2 = line([0, s*y(1)], [0, s*y(2)]);
set(h1, 'Color', c);
set(h2, 'Color', c);

function line_markers(x)
lim = get(gca,'YLim');
h = arrayfun(@(x) line([x x], lim), x);
set(h, 'Color', 'red');
