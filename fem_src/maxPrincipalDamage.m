classdef maxPrincipalDamage

    properties

        threshold; % value for damage threshold

    end

    methods

        function obj = maxPrincipalDamage(threshold)

            obj.threshold = threshold;

        end

        function D = damageCalc(obj,sig)

            global dim

%             if dim == 2
%
%                 stress = [sig(1) sig(3);...
%                           sig(3) sig(2)]; % 2D expanded stress tensor at node (2 x 2 MATRIX of FLOATs)
%
%             elseif dim == 3
%
%                 stress = [sig(1) sig(6) sig(5);...
%                           sig(6) sig(2) sig(4);...
%                           sig(5) sig(4) sig(3)]; % 3D expanded stress tensor at node (3 x 3 MATRIX of FLOATs)
%
%             end
%
%             [~,PrincStress] = eig(stress); % solve the eigenvalue problem to produce the principal stresses and their directions
%             PrincStress = diag(PrincStress); % convert to a vector
%             mP = max(PrincStress); % find the max principal stress

            if dim == 2

                solution1 = ((sig(1)+sig(2))/2) + sqrt(((sig(1)+sig(2))/2)^2 + sig(3)^2);
                solution2 = ((sig(1)+sig(2))/2) - sqrt(((sig(1)+sig(2))/2)^2 + sig(3)^2);

                mP = max([solution1 solution2]);

            elseif dim == 3

                stress = [sig(1) sig(6) sig(5);...
                          sig(6) sig(2) sig(4);...
                          sig(5) sig(4) sig(3)]; % 3D expanded stress tensor at node (3 x 3 MATRIX of FLOATs)
                [~,PrincStress] = eig(stress); % solve the eigenvalue problem to produce the principal stresses and their directions
                PrincStress = diag(PrincStress); % convert to a vector
                mP = max(PrincStress); % find the max principal stress

            end

            D = mP/obj.threshold; % calculate damage value

        end

    end

end
