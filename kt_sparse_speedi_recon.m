%
% Image reconstruction for the k-t sparse SPEEDI data acquisition using
% the joint partial spatial-temporal separablity
% (reference: Zhao et al., (2012) IEEE TMI 31:1809-1820)
%
% reconmode: reconstruction mode
%        = 1: continuous mode: all echoes are reconstructed together
%        = 2: periodic mode: echos are reconstructed in time blocks
% kdata: k-space raw data (5-D) [phase,freq,channel,echo,time lag]
%        for time-lock acquisition: time lag is the number of 
%          time delays (e.g., cardiac phases) relative to the
%          trigger event(e.g., cardiac gating).
%        for real-time acquisition: time lag is same as the
%          number of image time frames.
% kymode: k-y table file orgnization mode:
%        = 1: real-time mode
%        = 2: time-lock mode
% kyfile: ky index table file used for data acquisition
% reorglag: reorgnize the time lag order in kspace and ky file
%          =[]: no reorgnization
%          =[Nlag x 1]: new time lag order
% Np:   number of phase encodings in full sampling images
% Nnav: number of center k-space lines (navigators) for calculating
%       temporal subspace. These lines are sampled all time points,i.e.,
%       the center k-space are fully sampled
% Ny_recon: number of reconstructed pixels in phase encoding direction
%       if = []: same as Np
% Vt:  basis vectors for the temporal subspace
%      = []: automatic calculation using SVD on the navigator k-data
%      = [L x M x P]: predetermnined Vt
%        L: model order; M: number of time points
%        P: = 1 when reconmode = 1; = number of echoes when reconmode = 2      
% varargin{1}: order of ps model (L), default = 32
%              ONLY valid if Vt is not empty
% varargin{2}: number of echos per time block (default = 1)
%              valid only when reconmode = 2 only 
% recon_pssp: reconstructed images from the sparsely sampled data
%               if #channel = 1: complex image
%               if #channel > 1: channel-combined magnitude images 
% varargout{1}: indivdiual channel complex images
% ----------------------------------------------------------------------------------------------------
%
%                       Unversity of Illinois at Chicago
%
%  10/30/2020   Qingfei Luo     original
%
% ----------------------------------------------------------------------------------------------------

function [recon_pssp,varargout] = kt_sparse_speedi_recon(reconmode, kdata, kymode, kyfile, reorglag, Np, Nnav, Ny_recon, Vt, varargin)

% 
% parameters for the PSS recon algorithm
%
beta = 1e3;  % accuracy parameter (alpha in the paper)
mu = 2.5e-3; %2.5e-6; % regularization parameter ("lambda" in the paper)

%
% read the acquisition ky table file and reorgnize k-space data
%
[Nsam,Nf,Nchan,Necho,Nlag] = size(kdata); % Nsam is the number of sampled k-y lines in one image frame
Nfr  = Necho*Nlag;  % total time frames
npixecho = Np*Nf*Nchan;
if isempty(Ny_recon)
    Ny_recon = Np;
end
nfillky = Ny_recon-Np; % number of zero k-ys filled in the k-space
if Ny_recon > Nf
    nfillkx = Ny_recon-Nf;
else
    nfillkx = 0;
end
kytable = load(kyfile);
kytable = kytable(1:Nsam*Nfr);
if ~isempty(reorglag)
    kytable = reshape(kytable,Necho,Nlag,Nsam);
    kytable = kytable(:,reorglag,:);
    kytable = reshape(kytable,[],1);
    kdata = kdata(:,:,:,:,reorglag);
end
[kdata, Navdata, Mask, padk] = reorg_kdata(3, kdata, kymode, kytable, Np, Nnav);

% get temporal subspace
if isempty(Vt) % PCA (SVD)
    if isempty(varargin)
        r = 32;         % the order of PS model (L),i.e., use the first r temporal components as the principle components of temporal signal
    else
        r = varargin{1};
    end
%     [~, ~, Vt_r] = svds(Navdata, r); % PCA on navigator data
%     Vt_r = Vt_r';  % principle temporal signal components (temporal subspace basis)
else
    r = size(Vt, 1); % moder order
%     Vt_r = Vt;
end

% resampling data for time blocks
switch reconmode
    case 1
        ntblock = 1;  % recon time blocks
        nfrblock = Nfr;
        mask_block{1} = Mask;
        kdata_block{1} = kdata;
        navdata_block{1} = Navdata;
%         Vt_r_block = zeros(r,Nfr,1);
%         Vt_r_block(:,:,1) = Vt_r;
        
    case 2
        if ~isempty(varargin)
            if length(varargin)>1
                nechoblock = varargin{2}; % number of echos per time block
            end
        end
        ntblock = Necho-nechoblock+1;
        nfrblock = Nlag*nechoblock;
        mask_block = cell(ntblock,1);
        kdata_block = cell(ntblock,1);
        navdata_block = cell(ntblock,1);
%         Vt_r_block = zeros(r,nfrblock,ntblock);
        tmpmask = zeros(npixecho,nfrblock);
        tmppadk = zeros(Np,Nf,Nchan,nfrblock);
        tmpnavdata = zeros(Nnav*Nf*Nchan,nfrblock);
        for tbi=1:ntblock
            for j=0:nechoblock-1 
             tmpmask(:,j+1:nechoblock:end)= Mask(:,tbi+j:Necho:end);  
             tmppadk(:,:,:,j+1:nechoblock:end) = padk(:,:,:,tbi+j:Necho:end);
             tmpnavdata(:,j+1:nechoblock:end) = Navdata(:,tbi+j:Necho:end);
%              Vt_r_block(:,j+1:nechoblock:end,tbi) = Vt_r(:,tbi+j:Necho:end);
            end
            mask_block{tbi} = logical(tmpmask);
            navdata_block{tbi} = tmpnavdata;
            tmp1 = reshape(tmppadk,[],1);
            tmp2 = reshape(tmpmask,[],1);
            kdata_block{tbi} = tmp1(tmp2>0);
        end
end

% 
% PSS reconstruction
%
fprintf('Model order %d PS-Sparse reconstruction \n', r);
recon_pssp = zeros(npixecho,Nfr);
for tbi = 1:ntblock
    % determine the temporal subspace from the central k-space (navigator) data
    if isempty(Vt) % PCA (SVD)
        [~, ~, Vt_r] = svds(navdata_block{tbi}, r); % PCA on navigator data
        Vt_r = Vt_r';  % principle temporal signal components (temporal subspace basis)
    else
        Vt_r = squeeze(Vt(:,:,tbi));
    end
    % determine the spatial subspace
    Us_r0 = zeros(npixecho, r);
%     tmpvtr = squeeze(Vt_r_block(:,:,tbi));
    Us_r = ps_sparse_recon(kdata_block{tbi}, mask_block{tbi}, Us_r0, Vt_r, mu, beta, 'xf_sparse', Np, Nf, Nchan, nfrblock);
    recon_block = Us_r*Vt_r;
%     recon_block = Us_r*tmpvtr;
    switch reconmode
        case 1
            recon_pssp = recon_block;
        case 2
            fprintf('finished reconstructing time block %d\n',tbi);
            recon_pssp(:,tbi:Necho:end) = recon_block(:,1:nechoblock:end);
            if (ntblock<Necho)&&(tbi==ntblock)
                for j=1:Necho-ntblock
                    recon_pssp(:,ntblock+j:Necho:end) = recon_block(:,j+1:nechoblock:end);
                end
            end
    end
end

% output reconstructed images
if nfillky>0 % zero-padding k-space
    recon_pssp = reshape(recon_pssp,Np,Nf,Nchan,Necho,Nlag);
    fftpssp = fft2(recon_pssp); % transform pssp image to k-space
    fftpssp = padarray(fftpssp,nfillky/2,nfillkx/2);
    recon_pssp = ifft2(fftpssp);
end
if Nchan > 1
    recon_pssp = reshape(recon_pssp,Ny_recon,Nf,Nchan,Necho,Nlag);
    mag_pssp = zeros(Ny_recon,Nf,Necho,Nlag); % magnitude images
    %     mag_fs = squeeze(sum(image_xyt.*conj(image_xyt),3)); % SOS combine
    %     mag_pssp = squeeze(sum(recon_pssp.*conj(recon_pssp),3));
    for k = 1:Nlag
        tmpimg = mean(sqrt(recon_pssp(:,:,:,:,k) .* conj(recon_pssp(:,:,:,:,k))), 3); % chanel-averaged magnitude image
        mag_pssp(:,:,:,k) = flip(fftshift(tmpimg,2),1);  % shift zero-compoenent to image center
    end
    varargout{1} = flip(fftshift(recon_pssp,2), 1);
    recon_pssp = mag_pssp;
else
    recon_pssp = reshape(recon_pssp,Ny_recon,Nf,Necho,Nlag);
    recon_pssp = flip(fftshift(recon_pssp),1); 
    varargout{1} = recon_pssp;
end



