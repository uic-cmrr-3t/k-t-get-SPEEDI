%
%       Create k-y sampling masks from a k-y index table
%
%  To achieve the k-t undersampling, the central k-y lines (navigators) 
%  are fully sampled/acquired, and the outer  k-y lines are randomly 
%  undersampled at a time frame. The navigator signal is used 
%  to obain the basis of temporal space using SVD (i.e.,PCA),
%  and the spatial images are reconstructed using the JPSS method
%  (Zhao et al., (2012) IEEE TMI 31:1809-1820).
%
% kymode: k-y table orgnization mode:
%         = 1: real time mode
%         = 2: time-lock mode
%         see create_kytable.m for details
% kytable: k-y index table [Nx1] used for image acquisition
% Np: number of phase encodings (k_y) in the full sampling image
% Nsam: number of sampled k-y lines in the undersampling image
% Nnav: number of navigator lines
%       ignored if Nsam = Np (full sampling)
% Nfr: number of time frames
% kymask: k_y sampling masks [Np, Nt]
% navmask: navigator masks [Np, 1] (=[] for full sampling)
% varargout{1}: ky index table [Nsam, Nfr]
% ------------------------------------------------------------------------------------------
% 
% 7/22/2020     Qingfei Luo     original 
%
% ------------------------------------------------------------------------------------------

function [kymask, navmask, varargout] = create_samp_mask(kymode, kytable, Np, Nsam, Nnav, Nfr)

kymask = zeros(Np, Nfr); % k_y sampling mask
navmask = zeros(Np,1); % navigator sampling mask
nav_ind = (Np/2-Nnav/2+1:Np/2+Nnav/2)'; % # navigator k-y indices
navmask(nav_ind) = 1;
navmask = logical(navmask);
switch kymode
    case 1  % real-time mode
        kytable = reshape(kytable,Nsam,Nfr);
        
    case 2  % time-lock mode
        kytable = reshape(kytable, Nfr, Nsam);
        kytable = kytable';
end

for ti = 1:Nfr
    kymask(kytable(:,ti),ti) = 1;   
end
kymask = logical(kymask);
varargout{1} = kytable;

