classdef FEDEIM < KerMorObject
    
    properties
    end
    
    methods
        
        function err = getInterpolErrors(~, basis, pts, data, sigma)
            out = 3;
            if nargin < 5
                out = 6;
            end
            np = length(pts);
            err = zeros(out,np);
            pi = ProcessIndicator('Computing error for %d points on %dx%d data (fast)',np,false,np,size(data,1),size(data,2));
            for k=1:np
                P = sparse(pts(1:k),1:k,ones(k,1),size(basis,1),k);
                c = basis / (P'*basis);
                residual = data - c*(P'*data);
                err(1,k) = max(abs(residual(:)));
                err(2,k) = max(Norm.L2(residual));
                err(3,k) = mean(Norm.L2(residual));
                if nargin == 5
                    residual = sigma*residual;
                    err(4,k) = max(abs(residual(:)));
                    err(5,k) = max(Norm.L2(residual));
                    err(6,k) = mean(Norm.L2(residual));
                end
                pi.step;
            end
            pi.stop;
        end
        
        function plot(~, pm, tag, S, e, req_e)
%             ax = pm.nextPlot([tag '_sv'],['Singular value decay: ' tag],...
%                 'Size','Singular value');
%             semilogy(ax,1:size(S,2),diag(S));

            err = e(3,:);
            elow = find(err < 1e-8,1,'first');
            
            ax = pm.nextPlot([tag '_err'],...
                sprintf('Residual errors ''%s''\nerror < 1e^{-10} at m=%d, elements: %d', tag, elow, req_e(elow)),...
                'Size','Error');
            n = 1:size(e,2);
            ax2 = plotyy(ax,n,err,n,req_e,'semilogy','plot');
            linkaxes(ax2,'x');
            axis(ax2,'tight');
        end
        
        function [pts, v, elems] = getInterpolationPointskNN(~, u, k, idx_elements)
            n = size(u,1);
            M = size(u,2);

            % Init
            pts = zeros(1,M);
            v = pts;
            [v(1), pts(1)] = max(abs(u(:,1)));
            P = zeros(n,1);
            P(pts(1)) = 1;
            elems = find(idx_elements(:,pts(1)));
            i = 2;
            keff = k;
            while i < M+1
                imat = (P'*u(:,1:(i-1)));
                if keff/size(u,1) > .7 && rank(imat) < size(imat,1)
                    keff = keff-1;
                    i = i-1;
                else
                    c = imat \ (P'*u(:,i));
                    res = u(:,i) - u(:,1:(i-1))*c;
                    [val_cand, pts_cand] = sort(abs(res),'descend');
                    keff = k;
                end
                val_cand = val_cand(1:keff);
                pts_cand = pts_cand(1:keff);
                [aff_elem,~] = find(idx_elements(:,pts_cand));
                [~, ~, pos_aff] = intersect(elems,aff_elem,'stable');
                % Default: Make classic DEIM choice and use max abs value of residual
                sel = 1;
                if ~isempty(pos_aff)
                    % "k-neighbor" choice: If any point corresponding to the k largest
                    % residual values affects an element that is already used, take
                    % that point (and hence be a bit "less" optimal but use less
                    % elements)
                    sel = min(pos_aff);
                else
                    elems(end+1) = aff_elem(1);%#ok
                end
                pts(i) = pts_cand(sel);
                v(i) = val_cand(sel);
                P = sparse(pts(1:i),1:i,ones(i,1),n,i);
                i=i+1;
            end
        end
        
        function [pts, elems, basis_idx, maxerr, maxerr_res] = getInterpolationElements(this, u, idx_elems)
            sigma = double(idx_elems);
            num_elems = size(sigma,1);
            n = size(u,1);
            elems = zeros(1,num_elems);
            maxerr = zeros(1,num_elems+1);
            maxerr_res = maxerr;
            
            % Init
            resi_idx = 1:size(u,2);
            residual = u;
            pts = [];
            basis_idx = [];
            aborted = false;
            for elem_idx = 1:num_elems
                % Determine element error measure
                %     no = sqrt(sigma*(residual.^2));
                no = sigma*abs(residual);
                [maxerr(elem_idx), elems(elem_idx)] = max(max(no,[],2));
                maxerr_res(elem_idx) = max(abs(residual(:)));
                % Get corresponding points
                pts_elem = find(idx_elems(elems(elem_idx),:));
                P_elem = sparse(pts_elem,1:length(pts_elem),ones(length(pts_elem),1),n,length(pts_elem));
                
                % Do a local EIM point search on restricted set
                [pts_local, basis_local] = this.getInterpolationPointsMaxSearch(P_elem'*residual);
                % Bring the points in the order the've been chosen
                pts_elem = pts_elem(pts_local);
                % Add to global point set
                pts = [pts pts_elem];%#ok
                % Add points to basis (index is residual-dependent, so transfer to
                % global index)
                basis_idx = [basis_idx resi_idx(basis_local)];%#ok
                % Remove selected new basis from residual index set
                resi_idx(basis_local) = [];
                
                if isempty(resi_idx)
                    aborted = true;
                    break;
                end
                
                % Update residual
                basis = u(:,basis_idx);
                npts = length(pts);
                P = sparse(pts,1:npts,ones(npts,1),n,npts);
                c = (P'*basis) \ (P'*u(:,resi_idx));
                residual = u(:,resi_idx) - basis*c;
            end
            if aborted
                elems = elems(1:elem_idx);
            end
        end
        
        function [pts, basis_idx, v] = getInterpolationPointsMaxSearch(this, u)
            % Computes the interpolation indices according to the DEIM
            % algorithm
            %
            % Parameters:
            % u: matrix, columns are the orthonormal basis vectors computed via POD from the snapshot data
            % @type matrix<double>
            n = size(u,1);
            % Only search over            
            m = min(n,size(u,2));
            pts = zeros(1, n);
            v = pts;
            [val, idx]= max(abs(u));
            [v(1), j] = max(val);
            pts(1) = idx(j);
            P = zeros(n,1);
            P(pts(1)) = 1;
            basis_idx = zeros(1,m);
            resi_idx = 1:size(u,2);
            basis_idx(1) = j;
            resi_idx(j) = [];
            aborted = false;
            for k=2:m
                inv = (P'*u(:,basis_idx(1:k-1)));
                if rank(inv) < size(inv,2)
                    aborted = true;
                    k = k-1;%#ok
                    break;
                end
                c = inv \ (P'*u(:,resi_idx));
                resi = u(:,resi_idx) - u(:,basis_idx(1:k-1))*c;
                % Get absolute column maxima of residual (=~basis)
                [val, idx]= max(abs(resi),[],1);
                % Get maxima of all columns
                [v(k), j] = max(val);
                pts(k) = idx(j);
                % Save new basis vector
                basis_idx(k) = resi_idx(j);
                % Remove from remaining vectors to search
                resi_idx(j) = [];
                P = sparse(pts(1:k),1:k,ones(k,1),n,k);
                
                if isempty(resi_idx)
                    aborted = true;
                    break;
                end
            end
            if aborted
                pts = pts(1:k);
                basis_idx = basis_idx(1:k);
            end
        end
    end
    
end

% idx_elements = f.idx_vp_dof_unass_elems;
% sigma = double(idx_elems);
% num_elems = size(sigma,1);
% n = size(u,1);
% elems = zeros(1,num_elems);
% maxerr = zeros(1,num_elems+1);
% maxerr_res = maxerr;
% 
% % Init
% resi_idx = 1:size(u,2);
% resi_pos = true(1,size(u,2));
% residual = u;
% pts = [];
% basis_idx = [];
% aborted = false;
% npts = 0;
% for elem_idx = 1:num_elems
%     % Determine element error measure
%     %     no = sqrt(sigma*(residual.^2));
%     no = sigma*abs(residual);
%     
%     [maxerr(elem_idx), elems(elem_idx)] = max(max(no,[],2));
%     maxerr_res(elem_idx) = max(abs(residual(:)));
%     
%     % Get corresponding points
%     pts_elem = find(idx_elems(elems(elem_idx),:));
%     
%     % Do a local DEIM point search on restricted set
%     [v(end+1), pt_loc] = max(abs(u(pts_elem,npts+1)));
%     pts = [pts pts_elem(pt_loc)];%#ok
%     npts = length(pts);
%     P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     v_loc_max = v(end);
%     for i = npts + (2:length(pts_elem))
%         if v(end) < v_loc_max*sqrt(eps);
%             break;
%         end
%         c = (P'*u(:,1:(i-1))) \ (P'*u(:,i));
%         [v(end+1), pt_loc] = max(abs(u(pts_elem,i) - u(pts_elem,1:(i-1))*c));
%         pts = [pts pts_elem(pt_loc)];%#ok
%         npts = length(pts);
%         P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     end
%     resi_pos(pts) = false;
%     c = (P'*u(:,1:npts)) \ (P'*u(:,resi_pos));
%     residual = u(:,resi_pos) - u(:,1:npts)*c;
%     
%     %     [a,b,c] = svd(Psel'*hlp,'econ');
%     %     basis = [basis o.orthonormalize(u*c)];%#ok
%     %     pts = [pts newpts];%#ok
%     %     npts = length(pts);
%     %     P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     %     c = (P'*basis) \ (P'*u);
%     %     residual = u - basis*c;
%     %     hlp = u - basis*basis'*u;
% end
% if aborted
%     elems = elems(1:elem_idx);
% end
% err_FEDEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass);
% % err_FEDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass, s.f.Sigma);
% req_elems_FEDEIM = zeros(1,length(pts_FEDEIM));
% for k = 1:length(pts_FEDEIM)
%     req_elems_FEDEIM(k) = sum(sum(s.f.Sigma(:,pts_FEDEIM(1:k))*idx_elems(:,pts_FEDEIM(1:k))') ~= 0);
% end

