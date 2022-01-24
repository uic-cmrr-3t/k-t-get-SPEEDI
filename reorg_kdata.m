%
%  Reorgnize the raw k-space data according to the k-y index table
%
%  The raw data will be sorted and converted to 1-D data with
%  the order: phase, frequency, channel, echo, time lag.
%
%  If one ky line is acquired for multiple times at one frame,
%  then the mean value will be used to replace the raw data.
%
% acqmode: acquisition mode (see create_table.m for detals)
%        = 1: conventional spin/gradent echo
%        = 2: single-shot EPI
%        = 3: EPI-SPEEDI
% kdata: raw undersampled k-space data (5-D) [phase,freq,channel,echo,time lag]
% kymode: k-y table orgnization mode:
%         = 1: real time mode
%         = 2: time-lock mode
%         see create_kytable.m for details
% kytable: k-y index table [Nx1] used for image acquisition
% Np: number of phase encodings (k_y) in the full sampling image
% Nnav: number of navigator lines
%       Nnav is ignored for full sampling data
% newkdata: reorgnized k-space data (1-D) in order of PE,FE,Channel,Time
% navkdata: navigator k-space data [Nnav*Nf*Nchan,Nfr]
% kmask: k-space data sampling mask [Np*Nf*Nchan,Nfr]
% varargout{1}: raw k-space data with zero-paddings [Np, Nf, Nchan, Nfr]
% varargout{2}: ky sampling mask [Np, Nfr]
% varargout{3}: navigator sampling mask [Np, 1]
% varargout{4}: ky table [Nsam,Nfr]
% ------------------------------------------------------------------------------------------
% 
% 8/30/2020     Qingfei Luo     original 
%
% ------------------------------------------------------------------------------------------

function [newkdata, navkdata, kmask, varargout] = reorg_kdata(acqmode, kdata, kymode, kytable, Np, Nnav)

[Nsam,Nf,Nchan,Necho,Nlag] = size(kdata); % Nsam is the number of sampled k-y lines in one image frame
switch acqmode
    case 3 % EPI-SPEEDI
        Nfr  = Necho*Nlag;  % total number of image time frames
    otherwise
        Nfr = Nlag;
        if Necho>1
            fprintf('error: number of echos must be 1!\n');
            return;
        end
end

[kymask, navmask, kytable] = create_samp_mask(kymode, kytable, Np, Nsam, Nnav, Nfr);
kdata = reshape(kdata,Nsam,Nf,Nchan,Nfr);
kmask = repmat(kymask,[1,Nf]);
kmask = reshape(kmask,Np,Nfr,Nf);
kmask = permute(kmask,[1,3,2]); % k-space mask [Np,Nf,Nfr]
kmask = reshape(kmask,Np*Nf,Nfr);
kmask = repmat(kmask,[Nchan,1]); % k-space mask [Np*Nf*Nchan,Nfr]
newkdata = zeros(sum(kmask(:)),1);
padkdata = zeros(Np,Nf,Nchan,Nfr); % zero-padded k-space data
navkdata = zeros(Nnav,Nf,Nchan,Nfr);
navind = sort(find(navmask)); % navigator indice
if Nsam < Np  % sparse sampling
    kid = 0;
    for ti = 1:Nfr
        [sortky,oldid]=sort(kytable(:,ti));
        kdatati = kdata(oldid,:,:,ti); % sorted k-space data at this time frame
        [~,ia,~]=intersect(sortky,navind);
        navdatati = kdatati(ia,:,:,1);  % sorted navigator data at this time frame
        uqky = unique(sortky);  % sampling ky indices for this time frame
        nuqky = length(uqky);
        if nuqky<Nsam   % check if any ky are acquired for multiple times
            kdatati = zeros(nuqky,Nf,Nchan,1);
            for j=1:length(uqky) % fill the k-space and navigator data
                tmp = kdata(kytable(:,ti)==uqky(j),:,:,ti);
                kdatati(j,:,:,1)  = mean(tmp,1); % take the mean value if the ky has multiple acquisitions
                if ismember(uqky(j),navind)  % fill the navigator k-space dataset
                    navdatati(navind==uqky(j),:,:,1) = kdatati(j,:,:,1);
                end
            end
        end
        padkdata(uqky,:,:,ti) = kdatati;
        nptti = nuqky*Nf*Nchan;
        newkdata(kid+1:kid+nptti) = reshape(kdatati,nptti,1);
        navkdata(:,:,:,ti) = navdatati;
        kid = kid+nptti;
    end
    
else    % full sampling
    kid = 0;
    for ti = 1:Nfr
        [sortky,oldid]=sort(kytable(:,ti));
        kdatati = kdata(oldid,:,:,ti); % sorted k-space data at this time frame
        uqky = unique(sortky);  % sampling ky indices for this time frame
        nuqky = length(uqky);
        if nuqky<Nsam   % check if any ky are acquired for multiple times
            kdatati = zeros(nuqky,Nf,Nchan,1);
            for j=1:length(uqky) % fill the k-space and navigator data
                tmp = kdata(kytable(:,ti)==uqky(j),:,:,ti);
                kdatati(j,:,:,1)  = mean(tmp,1); % take the mean value if the ky has multiple acquisitions
            end
        end
        padkdata(uqky,:,:,ti) = kdatati;
        navkdata(:,:,:,ti) = padkdata(navind,:,:,ti);
        nptti = nuqky*Nf*Nchan;
        newkdata(kid+1:kid+nptti) = reshape(kdatati,nptti,1);
        kid = kid+nptti;
    end
end
navkdata = reshape(navkdata,Nnav*Nf*Nchan,Nfr);
varargout{1} = padkdata;
varargout{2} = kymask;
varargout{3} = navmask;
varargout{4} = kytable;




