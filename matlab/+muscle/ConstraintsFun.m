classdef ConstraintsFun < dscomponents.ACoreFun
    %CONSTRAINTSFUN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function this = ConstraintsFun(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.update;
        end
        
        function update(this)
            [~,JP] = this.System.f.computeSparsityPattern;
            this.JSparsityPattern = JP;
            this.fDim = this.System.NumAlgebraicDofs;
            this.xDim = this.System.NumAlgebraicDofs;
        end
        
        function J = getStateJacobianImpl(this, ~, ~)
            f = this.System.f;
            J = f.curJGC;
            if isempty(J)
                error('Invalid use. This is a efficiency hack, not a standalone function.');
            else
%                 f.curJGC = [];
            end
        end
        
        function evaluateCoreFun(~)
            error('Do not call. evaluate is implemented directly.');
        end
        
        function gc = evaluate(this, ~, ~)
            f = this.System.f;
            gc = f.curGC;
            if isempty(gc)
                error('Invalid use. This is a efficiency hack, not a standalone function.');
            else
%                 f.curGC = [];
            end
        end
        
        function projected = project(this, V, W)
            projected = this.clone;
            projected = project@dscomponents.ACoreFun(this, V, W, projected);
        end
        
        function copy = clone(this)
            % Create new instance
            copy = muscle.ConstraintsFun(this.System);
            copy = clone@dscomponents.ACoreFun(this, copy);
        end
    end
    
end

