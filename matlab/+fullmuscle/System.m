classdef System < muscle.System;
% System: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties
        
    end
    
    properties(SetAccess=private)
        num_motoneuron_dof;
        num_sarco_dof;
        off_moto;
        off_sarco;
        
        % Membrane capacitance. Different values for different fibre types, 
        % due to different action potential propagation speeds
        % C_m is computed parameter-dependent.
        % These constants are for both slow and fast muscle models and are also used in the
        % first entry of the sarcomere constants computed in
        % models.muscle.FibreDynamics.initSarcoConst @type double
        C_m_slow = 0.58;
        C_m_fast = 1;
    end
    
    properties(Access=private)
        nfibres;
    end
    
    methods
        function this = System(model)
            this = this@muscle.System(model);
            this.nfibres = length(model.FibreTypes);
            this.f = fullmuscle.Dynamics(this);
        end
        
        function plot(this, t, y, varargin)
            plot@muscle.System(this, t, y(1:this.num_uvp_dof,:), varargin{:});
        end
        
%         function configUpdated(this)
%             configUpdated@muscle.System(this);
%              
%         end
    end
    
    methods(Access=protected)
        function updateDofNums(this, mc)
            updateDofNums@muscle.System(this, mc);
            
            this.num_motoneuron_dof = 6*this.nfibres;
            % Motoneurons are beginning after mechanics
            this.off_moto = this.num_uvp_dof; 

            % Sarcomeres are beginning after motoneurons
            this.num_sarco_dof = 56*this.nfibres;
            this.off_sarco = this.off_moto + this.num_motoneuron_dof;
        end
        
        function x0 = assembleX0(this)
            x0 = zeros(this.num_uvp_dof ...
                + this.num_motoneuron_dof + this.num_sarco_dof,1);
            % Get muscle x0
            x0(1:this.num_uvp_dof) = assembleX0@muscle.System(this);
            
            % Load x0 coefficients for moto/sarco system from file
            mc = metaclass(this);
            s = load(fullfile(fileparts(which(mc.Name)),'x0coeff.mat'));
            x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            
            ft = this.Model.FibreTypes;
            for k=1:this.nfibres
                x0ms = x0_motorunit.evaluate(ft(k));
                % add moto
                x0(this.off_moto + 6*(k-1) + (1:6)) = x0ms(1:6);
                % add sarco
                x0(this.off_sarco + 56*(k-1) + (1:56)) = x0ms(7:end);
            end
        end
        
        function Daff = assembleDampingMatrix(this)
            Daff_mech = assembleDampingMatrix@muscle.System(this);
            
            Daff = dscomponents.AffLinCoreFun(this);
            extra = (6+56)*this.nfibres;
            D = blkdiag(Daff_mech.getMatrix(1),sparse(extra,extra));
            Daff.addMatrix('mu(1)',D);
        end
        
        function M = assembleMassMatrix(this)
            M = assembleMassMatrix@muscle.System(this);
            extra = (6+56)*this.nfibres;
            M = blkdiag(M,speye(extra));
        end
    end
    
end