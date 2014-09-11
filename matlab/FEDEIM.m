classdef FEDEIM < KerMorObject
    
    properties
    end
    
    methods
        
        function err = getInterpolErrors(this, basis, pts, data, sigma)
            if nargin < 5
                sigma = 1;
            end
            np = length(pts);
            err = zeros(2,np);
            pi = ProcessIndicator('Computing error for %d points on %dx%d data',np,false,np,size(data,1),size(data,2));
            for k=1:np
                P = sparse(pts(1:k),1:k,ones(k,1),size(basis,1),k);
                c = (P'*basis) \ (P'*data);
                residual = sigma*(data - basis*c);
                err(1,k) = max(abs(residual(:)));
                err(2,k) = max(Norm.L2(residual));
                pi.step;
            end
            pi.stop;
        end
        
        function plot(this, pm, tag, S, e)
            ax = pm.nextPlot([tag '_sv'],['Singular value decay: ' tag],...
                'Size','Singular value');
            semilogy(ax,1:size(S,2),diag(S));
            ax = pm.nextPlot([tag '_err'],['Residual errors: ' tag],...
                'Size','Error');
            n = 1:size(e,2);
            semilogy(ax,n,e(1,:),'r',n,e(2,:));
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
                c = (P'*u(:,basis_idx(1:k-1))) \ (P'*u(:,resi_idx));
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
            if aborted > 0
                pts = pts(1:k);
                basis_idx = basis_idx(1:k);
            end
        end
    end
    
end

