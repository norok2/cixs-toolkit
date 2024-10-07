% set experment specific parameters for data analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% raw data directory
data_dir = [pwd '/../data/'];
plot_dir = [pwd '/../figures/'];
import_dir = [pwd '/../import/'];
export_dir = [pwd '/../export/'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% physical constants
% Planck constant h in J s
[h_P, UM_h_P, D_h_P] = physical_constant('PLANCK_CONSTANT');
% Speed of Light (Photon) in Vacuum c in m / s
[c_p, UM_c_p, D_c_p] = physical_constant('SPEED_OF_LIGHT_IN_VACUUM');
% Boltzmann constant k_B in J / K
[k_B, UM_k_B, D_k_B] = physical_constant('BOLTZMANN_CONSTANT');
% Elementary charge e in C
[q_e, UM_q_e, D_q_e] = physical_constant('ELEMENTARY_CHARGE');
% classical electron radius (r_e = e^2 / m_e / c^2) in m
[r_e, UM_r_e, D_r_e] = physical_constant('CLASSICAL_ELECTRON_RADIUS');
% mass of the electron in in kg
[m_e, UM_m_e, D_m_e] = physical_constant('ELECTRON_MASS');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimental parameters
% estimated experimental room temperature in K
expt.T = 20 + 273.15;
% base energy of monochromatized beam (corrections: #1: logbook note, #2: elastic centered to zero)
expt.E_in = 7.9070e3 + 0.375;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% beam parameters
% polarization fraction (horizontal / vertical) or (sigma / pi)
beam.polarization_fraction = 0.99;
beam.polarization = [beam.polarization_fraction, (1 - beam.polarization_fraction)];
% beam angular divergence in rad
beam.divergence = 3e-6;

% fundamental direction of experimental setup
beam.direction = [1, 0, 0];
beam.normal = [0, 0, 1]; % within beam and bragg-reflection plane
beam.binormal = [0, 1, 0]; % orthogonal to beam and bragg-reflection plane


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% monochromator (Silicon, Si using [111] reflection)
monochr.name = 'Silicon';
monochr.symbol = 'Si';
monochr.type = 'crystal';
monochr.element = [14];
% lattice parameters of Si in m
monochr.cell_size = [0.5431020504e-9; 0.5431020504e-9; 0.5431020504e-9;]; % (error: 0.0000000089e-9)
monochr.cell_angle = deg2rad([90, 90, 90]);
% position of atoms in crystal
monochr.unit_cell = [ % x, y, z, occupancy, Z_element
	(0/4), (0/4), (0/4), 1, 14; 
	(2/4), (2/4), (0/4), 1, 14;
	(2/4), (0/4), (2/4), 1, 14;
	(0/4), (2/4), (2/4), 1, 14;
	(1/4), (1/4), (1/4), 1, 14;
	(3/4), (3/4), (1/4), 1, 14;
	(3/4), (1/4), (3/4), 1, 14;
	(1/4), (3/4), (3/4), 1, 14;
	];
monochr.cell_offset = [(1/8), (1/8), (1/8)]; % shift the origin to inversion center
% density of Si in kg/m^3
monochr.density = 2.33e3;
% asymmetry angle of the sample
monochr.asymmetry = deg2rad(0.00+0.00);
% linear size in m
monochr.thickness = 0.1;
% temperature in K
monochr.T = 77;
% pressure in Pa
monochr.P = 1e5;
% used reflection
monochr.hkl = [1, 1, 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sample parameters (Silicon, Si)
sample.name = 'Silicon';
sample.symbol = 'Si';
sample.type = 'crystal';
% lattice parameters of Si in m
sample.cell_size = [0.5431020504e-9; 0.5431020504e-9; 0.5431020504e-9;]; % (error: 0.0000000089e-9)
sample.cell_angle = deg2rad([90, 90, 90]);
% position of atoms in crystal
sample.unit_cell = [ % x, y, z, occupancy, elem_id
	(0/4), (0/4), (0/4), 1, 14; 
	(2/4), (2/4), (0/4), 1, 14;
	(2/4), (0/4), (2/4), 1, 14;
	(0/4), (2/4), (2/4), 1, 14;
	(1/4), (1/4), (1/4), 1, 14;
	(3/4), (3/4), (1/4), 1, 14;
	(3/4), (1/4), (3/4), 1, 14;
	(1/4), (3/4), (3/4), 1, 14;
	];
sample.cell_offset = [(1/8), (1/8), (1/8)]; % shift the origin to inversion center
% density in kg/m^3
sample.density = 2.33e3;
% asymmetry angle of the sample
sample.asymmetry = deg2rad(5.00+0.00);
% sample linear size in m
sample.thickness = 50e-3;
% temperature of the sample in K
sample.T = 300;
% pressure of the sample in Pa
sample.P = 1e5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% energy dependent parameters (and related) for base energy value 
% Miller indexes of selected Bragg reflection
hkl = monochr.hkl;
E = expt.E_in;
crystal = monochr;
% Bragg angle in rad
th_B = bragg_angle(crystal, E, hkl);
% structure form factor in forward direction
[F_0, chi_0] = structure_form_factor(crystal, E, [0 0 0]);
% structure form factor in [hkl] direction
[F_h, chi_h] = structure_form_factor(crystal, E, hkl, +1, +1, -1);
% structure form factor in -[hkl] direction
[F_hb, chi_hb] = structure_form_factor(crystal, E, -hkl);

%% extra parameters
param.dE.min = 0; % eV
param.dE.max = 100; % eV
param.dE.step = 0.5;
% WARNING: note that this influences many things!!!
param.dth.min = -1e-2; % rad
param.dth.max = 1e-2; % rad
param.dth.sampling_ratio = 1;

param.plot.dE.min = -2; % eV
param.plot.dE.max = 60; % eV
param.plot.dth.min = -0.6e-4; % rad
param.plot.dth.max = 1.6e-4; % rad
param.plot.n_seq = 8;
param.plot.mesh = 100;

param.E.bound = 25; % eV
param.E.density = 4000;
param.th.bound = 1e-4; % rad
param.th.density = 4000;
param.precise_convolution = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experiment data information
expt_num = 8; % to keep updated with below

% build expt template
expt.param = param;
expt.beam = beam;
expt.monochr = monochr;
expt.sample = sample;
% supported: [selected_average/[2|3|4] | weighted_average/[2|3|4] | full_average/[2|3|4]]
expt.analysis = 'selected_average/2';
expt.name = '';
expt.signal.file = '';
expt.signal.scan = [];
expt.background.file = '';
expt.background.file = [];
expt.detector = '';
expt.group = 0;
expt = repmat(expt, 1, expt_num);


% see ID16 Hi-Res Logbook 108 R, p. 62 (14-20/02/2007)
% hkl456_2: [2,214; 447,659; 215,349; 350,445]; % [run 1, run 2, energy tails, background]
% hkl113: [2,274; 492,764; 275,394; 395,490]; % [run 1, run 2, energy tails, background]
% hkl327: [2,274; 492,764; 275,394; 395,490]; % [run 1, run 2, energy tails, background]
% hkl327_2: [2,274; 275,394]; % [run 1, energy tails]

expt(1).name = 'hkl456_1';
expt(1).signal.file = [data_dir, 'hkl456_2'];
expt(1).signal.scan = [2:214, 215:349];
expt(1).background.file = [data_dir, 'hkl456_2'];
expt(1).background.scan = [350:445];
expt(1).sample.hkl = [1, 1, 1];
expt(1).sample.phi = deg2rad(0.24);
expt(1).sample.asymmetry = geometric_asymmetry(expt(1).sample.asymmetry, expt(1).sample.phi);
expt(1).sample.q = -0.1 .* [5, 6, 4];
expt(1).sample.n_terms = 2;
expt(1).detector.distance = 1.041;
expt(1).detector.height = 0.240 + 0.028;
expt(1).detector.th = atan(expt(1).detector.height / expt(1).detector.distance);
expt(1).detector.phi = deg2rad(1.55);

expt(2).name = 'hkl456_2';
expt(2).signal.file = [data_dir, 'hkl456_2'];
expt(2).signal.scan = [447:659, 215:349];
expt(2).background.file = [data_dir, 'hkl456_2'];
expt(2).background.scan = [350:445];
expt(2).sample.hkl = [1, 1, 1];
expt(2).sample.phi = deg2rad(0.24);
expt(2).sample.asymmetry = geometric_asymmetry(expt(2).sample.asymmetry, expt(2).sample.phi);
expt(2).sample.q = -0.1 .* [5, 6, 4];
expt(2).sample.n_terms = 2;
expt(2).detector.distance = 1.041;
expt(2).detector.height = 0.240 + 0.028;
expt(2).detector.th = atan(expt(2).detector.height / expt(2).detector.distance);
expt(2).detector.phi = deg2rad(1.55);

expt(3).name = 'hkl113_1';
expt(3).signal.file = [data_dir, 'hkl113'];
expt(3).signal.scan = [2:274, 275:394];
expt(3).background.file = [data_dir, 'hkl113'];
expt(3).background.scan = [395:490];
expt(3).sample.hkl = [1, 1, 1];
expt(3).sample.phi = deg2rad(-2.00);
expt(3).sample.asymmetry = geometric_asymmetry(expt(3).sample.asymmetry, expt(3).sample.phi);
expt(3).sample.q = -0.5 .* [-1, 3, 1];
expt(3).sample.n_terms = 2;
expt(3).detector.distance = 1.041;
expt(3).detector.height = 0.240 + 0.005;
expt(3).detector.th = atan(expt(3).detector.height / expt(3).detector.distance);
expt(3).detector.phi = deg2rad(24.49);

expt(4).name = 'hkl113_2';
expt(4).signal.file = [data_dir, 'hkl113'];
expt(4).signal.scan = [492:764, 275:394];
expt(4).background.file = [data_dir, 'hkl113'];
expt(4).background.scan = [395:490];
expt(4).sample.hkl = [1, 1, 1];
expt(4).sample.phi = deg2rad(-2.00);
expt(4).sample.asymmetry = geometric_asymmetry(expt(4).sample.asymmetry, expt(4).sample.phi);
expt(4).sample.q = -0.5 .* [-1, 3, 1];
expt(4).sample.n_terms = 2;
expt(4).detector.distance = 1.041;
expt(4).detector.height = 0.240 + 0.005;
expt(4).detector.th = atan(expt(4).detector.height / expt(4).detector.distance);
expt(4).detector.phi = deg2rad(24.49);

expt(5).name = 'hkl327_1';
expt(5).signal.file = [data_dir, 'hkl327'];
expt(5).signal.scan = [2:274, 275:394];
expt(5).background.file = [data_dir, 'hkl327'];
expt(5).background.scan = [395:490];
expt(5).sample.hkl = [1, 1, 1];
expt(5).sample.phi = deg2rad(+1.75);
expt(5).sample.asymmetry = geometric_asymmetry(expt(5).sample.asymmetry, expt(5).sample.phi);
expt(5).sample.q = -0.25 .* [-3, 7, 2];
expt(5).sample.n_terms = 2;
expt(5).detector.distance = 1.041;
expt(5).detector.height = 0.240 + 0.001;
expt(5).detector.th = atan(expt(5).detector.height / expt(5).detector.distance);
expt(5).detector.phi = deg2rad(30.76);

expt(6).name = 'hkl327_2';
expt(6).signal.file = [data_dir, 'hkl327'];
expt(6).signal.scan = [492:764, 275:394];
expt(6).background.file = [data_dir, 'hkl327'];
expt(6).background.scan = [395:490];
expt(6).sample.hkl = [1, 1, 1];
expt(6).sample.phi = deg2rad(+1.75);
expt(6).sample.asymmetry = geometric_asymmetry(expt(6).sample.asymmetry, expt(6).sample.phi);
expt(6).sample.q = -0.25 .* [-3, 7, 2];
expt(6).sample.n_terms = 2;
expt(6).detector.distance = 1.041;
expt(6).detector.height = 0.240 + 0.001;
expt(6).detector.th = atan(expt(6).detector.height / expt(6).detector.distance);
expt(6).detector.phi = deg2rad(30.76);

%  expt(7).name = 'hkl327_3';
%  expt(7).signal.file = [data_dir, 'hkl327_2'];
%  expt(7).signal.scan = [2:274, 275:394];
%  expt(7).background.file = [data_dir, 'hkl327'];
%  expt(7).background.scan = [395:490];
%  expt(7).sample.hkl = [1, 1, 1];
%  expt(7).sample.phi = deg2rad(+1.75);
%  expt(7).sample.asymmetry = geometric_asymmetry(expt(7).sample.asymmetry, expt(7).sample.phi);
%  expt(7).sample.q = -0.25 .* [-3, 7, 2];
%  expt(7).detector.distance = 1.041;
%  expt(7).detector.height = 0.240 + 0.001;
%  expt(7).detector.th = atan(expt(7).detector.height / expt(7).detector.distance);
%  expt(7).detector.phi = deg2rad(30.76);

expt(7).name = 'hkl111_1';
expt(7).signal.file = [data_dir, 'hkl111'];
expt(7).signal.scan = [2:163, 188:298, 299:418];
expt(7).background.file = [data_dir, 'hkl111'];
expt(7).background.scan = [419:514];
expt(7).sample.hkl = [2, 2, 0];
expt(7).sample.phi = deg2rad(+0.00);
expt(7).sample.asymmetry = geometric_asymmetry(expt(7).sample.asymmetry, expt(7).sample.phi);
expt(7).sample.q = 1.32 .* [1, 1, 1];
expt(7).sample.n_terms = 3;
expt(7).detector.distance = 1.004;
expt(7).detector.height = 0.527;
expt(7).detector.th = atan(expt(7).detector.height / expt(7).detector.distance);
expt(7).detector.phi = deg2rad(25.96);


expt(8).name = 'hkl111_2';
expt(8).signal.file = [data_dir, 'hkl111'];
expt(8).signal.scan = [516:653, 140:163, 170:184, 203:298, 299:418];
expt(8).background.file = [data_dir, 'hkl111'];
expt(8).background.scan = [419:514];
expt(8).sample.hkl = [2, 2, 0];
expt(8).sample.phi = deg2rad(+0.00);
expt(8).sample.asymmetry = geometric_asymmetry(expt(8).sample.asymmetry, expt(8).sample.phi);
expt(8).sample.q = 1.32 .* [1, 1, 1];
expt(8).sample.n_terms = 3;
expt(8).detector.distance = 1.004;
expt(8).detector.height = 0.527;
expt(8).detector.th = atan(expt(8).detector.height / expt(8).detector.distance);
expt(8).detector.phi = deg2rad(25.96);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% configurations to be made after experiments definition
groups.n = 0;
for i=1:expt_num
	if (expt(i).group == 0)
		groups.n = groups.n + 1;
		for j=i:expt_num
			if (strcmp(expt(i).name(4:6), expt(j).name(4:6)))
				expt(j).group = groups.n;
			end
		end
	end
end

groups.dim = zeros(1, groups.n);
for i=1:expt_num
	groups.dim(expt(i).group) = groups.dim(expt(i).group) + 1;
end