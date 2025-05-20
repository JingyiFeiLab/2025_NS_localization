%%
% Outputs description:
% Cell_Summary: quantities at the single-cell level
% Cell_Summary_f (filted by removing outliers)
% spec_prop: quantities at the single-speckle level
% spec_prop_f (filted by removing outliers)
% bg: background information
% NSmorphology: quantities of nuclear speckles morphology
% NSmorphology_cell: quantities of nuclear speckles morphology at the
% single-cell level
load([cellname filesep consname 'cell_intensity.mat'])
load([nucleusname filesep consname 'nu_intensity.mat'])
load([spec_trans_name filesep consname 'spec_intensity.mat'])
load([cellname filesep consname 'cam_bg.mat'])
load([transname filesep consname '_trans.mat'])
cell_num_all=size(spec_intensity,1);

%% background estimation
bg(:,1)=bg(:,1)*remove_bg_R;
if remove_bg_G == 1
    bg(:,2)=bg(:,2)*remove_bg_G;
end
if remove_bg_G == 0
    bg(:,2)=bg(:,2)*remove_bg_G;
end
if remove_bg_G == 2
    bg(:,2)=cam_bg_G*ones(3,1);
end
bg(:,3)=bg(:,3)*remove_bg_B;

%% Create the summary table
Cell_Summary = table(zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), ...
    zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), ...
    zeros(cell_num_all, 3), zeros(cell_num_all, 3), ...
    zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), zeros(cell_num_all, 3), ...
    zeros(cell_num_all, 1), zeros(cell_num_all, 1), ...
    zeros(cell_num_all, 1), zeros(cell_num_all, 1), zeros(cell_num_all, 1), zeros(cell_num_all, 1), ...
    zeros(cell_num_all, 1), zeros(cell_num_all, 1), zeros(cell_num_all, 1), ...
    zeros(cell_num_all, 1), zeros(cell_num_all, 1), zeros(cell_num_all, 1), zeros(cell_num_all, 1), ...
    'VariableNames', {'I_cell', 'c_cell', 'I_cyto', 'c_cyto', ...
    'I_nu', 'c_nu', 'I_nu_r','c_nu_r',...
    'I_spec', 'c_spec', ...
    'I_nucleoplasm','c_nucleoplasm', 'I_nucleoplasm_r','c_nucleoplasm_r',...
    'id_cell','spec_vol_frac', ...
    'Area_cell','Area_nu','Area_nu_r','Area_spec', ...
    'A_nucleoplasm','A_nucleoplasm_r', 'NS_percell', ...
    'P_R','P_G', 'cyto_R','cyto_G'});

%% Compute total or mean intensities
cell_id=spec_intensity(:,5);
for i=1:cell_num_all
    id=cell_id(i);
    A_cell=cell_intensity(id,1);
    A_nu=nu_intensity(id,1);
    A_nu_r=nu_intensity_r(id,1);
    A_cyto=A_cell-A_nu;
    A_spec=spec_intensity(i,1);
    A_nucleoplasm=A_nu-A_spec;
    A_nucleoplasm_r=A_nu_r-A_spec;

    for j=1:3
        Cell_Summary.I_cell(i,j)=cell_intensity(id,j+1)-A_cyto*bg(3,j)-A_nu*bg(2,j);
        Cell_Summary.c_cell(i,j)=Cell_Summary.I_cell(i,j)/A_cell;

        Cell_Summary.I_cyto(i,j)=cell_intensity(id,j+1)-nu_intensity(id,j+1)-A_cyto*bg(3,j);
        Cell_Summary.c_cyto(i,j)=Cell_Summary.I_cyto(i,j)/A_cyto;

        Cell_Summary.I_nu(i,j)=nu_intensity(id,j+1)-A_nu*bg(2,j);
        Cell_Summary.c_nu(i,j)=Cell_Summary.I_nu(i,j)/A_nu;

        Cell_Summary.I_nu_r(i,j)=nu_intensity_r(id,j+1)-A_nu_r*bg(2,j);
        Cell_Summary.c_nu_r(i,j)=Cell_Summary.I_nu_r(i,j)/A_nu_r;

        Cell_Summary.I_spec(i,j)=spec_intensity(i,j+1)-A_spec*bg(2,j);
        Cell_Summary.c_spec(i,j)=Cell_Summary.I_spec(i,j)/A_spec;

        Cell_Summary.I_nucleoplasm(i,j)=Cell_Summary.I_nu(i,j)-Cell_Summary.I_spec(i,j);
        Cell_Summary.c_nucleoplasm(i,j)=Cell_Summary.I_nucleoplasm(i,j)/A_nucleoplasm;

        Cell_Summary.I_nucleoplasm_r(i,j)=Cell_Summary.I_nu_r(i,j)-Cell_Summary.I_spec(i,j);
        Cell_Summary.c_nucleoplasm_r(i,j)=Cell_Summary.I_nucleoplasm_r(i,j)/A_nucleoplasm_r;
    end

    Cell_Summary.id_cell(i,1)=id;
    Cell_Summary.spec_vol_frac(i,1)=A_spec/A_nu;
    Cell_Summary.NS_percell(i,1)=spec_intensity(i,6);
    Cell_Summary.Area_cell(i,1)=A_cell;
    Cell_Summary.Area_nu(i,1)=A_nu;
    Cell_Summary.Area_nu_r(i,1)=A_nu_r;
    Cell_Summary.Area_spec(i,1)=A_spec;
    Cell_Summary.A_nucleoplasm(i,1)=A_nucleoplasm;
    Cell_Summary.A_nucleoplasm_r(i,1)=A_nucleoplasm_r;
    Cell_Summary.P_R(i,1)=Cell_Summary.c_spec(i,1)/Cell_Summary.c_nucleoplasm_r(i,1);
    Cell_Summary.P_G(i,1)=Cell_Summary.c_spec(i,2)/Cell_Summary.c_nucleoplasm_r(i,2);
    Cell_Summary.cyto_R(i,1)=Cell_Summary.I_cyto(i,1)/Cell_Summary.I_cell(i,1);
    Cell_Summary.cyto_G(i,1)=Cell_Summary.I_cyto(i,2)/Cell_Summary.I_cell(i,2);
end

% Note: Cell_Summary is a table, while spec_prop is a structure array
% Use rmoutliers function to remove outliers
PR_raw=Cell_Summary.P_R;
[~, outlierIdx] = rmoutliers(PR_raw);
Cell_Summary_f = Cell_Summary(~outlierIdx, :);
spec_prop_f = spec_prop(~outlierIdx);
% Generate quantities in speckle morphology
all_v = arrayfun(@(x) transpose(x.singlecell.v), spec_prop_f, 'UniformOutput', false);
all_cir = arrayfun(@(x) transpose(x.singlecell.Circularity), spec_prop_f, 'UniformOutput', false);
all_ecc = arrayfun(@(x) transpose(x.singlecell.Eccentricity), spec_prop_f, 'UniformOutput', false);
all_v_T = transpose(all_v);
all_cir_T = transpose(all_cir);
all_ecc_T = transpose(all_ecc);
NSmorphology.Area = vertcat(all_v_T{:});
NSmorphology.Circularity = vertcat(all_cir_T{:});
NSmorphology.Eccentricity = vertcat(all_ecc_T{:});
% Generate quantities in speckle morphology at the single-cell level
NSmorphology_cell.Area = transpose(arrayfun(@(x) mean(x.singlecell.v), spec_prop_f));
NSmorphology_cell.Circularity = transpose(arrayfun(@(x) mean(x.singlecell.Circularity), spec_prop_f));
NSmorphology_cell.Eccentricity = transpose(arrayfun(@(x) mean(x.singlecell.Eccentricity), spec_prop_f));
% 
save([summary_folder filesep consname '_specP.mat'],'Cell_Summary','Cell_Summary_f','spec_prop','spec_prop_f','NSmorphology','NSmorphology_cell','bg')