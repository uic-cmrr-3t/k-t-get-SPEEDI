%
% do cross correlation b/w sparse and full sampling k-space phase lines
% 
% fullk: 3-D full sampling k-space data [Npe, Nfe, Nch]
% sparsek: 3-D sparse sampling k-space data [Ns, Nfe, Nch]
% cf: correlation coefficients b/w sparse and full sampling data  [Ns x Npe]
% varargout{1}: phase line indices in full sampling data that have the max
%               correlation with the sparse sampling phase lines
% varargout{2}: reorgnized and padded sparse sampling data [Npe, Nfe, Nch]
% ---------------------------------------------------------------------------------

function [cf, varargout] = kspace_xcorr(fullk, sparsek)

[Np, Nf, Nc] = size(fullk);
[Ns, ~, ~] = size(sparsek);
cf = zeros(Ns,Np); % correlation coefficients b/w sparse and full sampling data
maxpid = zeros(Ns,1);  % the phase line indices of full sampling data with max correlation
reorgsk = zeros(Np,Nf,Nc);
for sid = 1:Ns
    tmpcc = zeros(Nc,Np);
    sigs = squeeze(sparsek(sid,:,:));
    for pid = 1:Np
        for j=1:Nc
            tmp = corrcoef(abs(sigs(:,j)),abs(squeeze(fullk(pid,:,j))));
            tmpcc(j,pid) = tmp(1,2);
        end
    end
    tmpcc(isnan(tmpcc)) = 0;
    tmpm = mean(tmpcc,1);
    [~,tmpid] = max(tmpm);
    maxpid(sid) = tmpid;
    cf(sid,:) = tmpm;
    reorgsk(tmpid,:,:) = sigs;
end
varargout{1} = maxpid;
varargout{2} = reorgsk;