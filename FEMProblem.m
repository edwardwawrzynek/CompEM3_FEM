% 
classdef FEMProblem
    properties
        ep_r_center % pde parameter (alpha)
        beta  % pde parameter
        gamma % neumann bc term
        kappa_n % neumann bc term
        f % forcing
        geo_a % inner radius of domain
        geo_b % outer radius of domain

        model % pde model object
        mesh % pde mesh object
        K % system matrix
        b % right hand side
    end
    methods
        function obj = FEMProblem(ep_r_center, beta, gamma, kappa_n, f, geo_a, geo_b)
            obj.ep_r_center = ep_r_center;
            obj.beta = beta;
            obj.gamma = gamma;
            obj.kappa_n = kappa_n;
            obj.f = f;
            obj.geo_a = geo_a;
            obj.geo_b = geo_b;
        end
        
        % create the geometry and mesh the problem
        function obj = createMesh(obj, h)
            % geometry codes for pde toolbox
            circle_code=1;
            circle_footer=[];

            % create model
            obj.model = createpde();
            % outer circle (boundary)
            B = [circle_code, 0, 0, obj.geo_b, circle_footer]';
            % inner circle (scatterer)
            A = [circle_code, 0, 0, obj.geo_a, circle_footer]';
            set_names = ['B'; 'A']';
            % TODO: add parameter to mesh and tag (rather than exclude) scatterer (for dielectric case)
            set_formula = 'B - A';

            % build geometry
            gd = [B A];
            g = decsg(gd, set_formula, set_names);
            geometryFromEdges(obj.model, g);
            
            % mesh
            obj.mesh = generateMesh(obj.model, Hmax=h*2, Hmin=h, Hgrad=2, GeometricOrder="linear");
        end
        
        % plot the mesh
        % if doLabels is true the elements and nodes are labeled by index
        function plotMesh(obj, doLabels)
            if nargin < 2
                doLabels = true;
            end
            % nodes
            pos = obj.mesh.Nodes;
            % connectivity matrix
            n = obj.mesh.Elements;
            % number of nodes and triangles
            Nnod = length(pos);
            Ntri = length(n);
            
            figure();
            pdemesh(obj.model);
            if doLabels
                for i=1:Nnod
                    text(pos(1,i), pos(2,i), num2str(i), "HorizontalAlignment", "center", "VerticalAlignment","middle", "BackgroundColor", [0.8,0.8,0.8,1]);
                end
                for i=1:Ntri
                    text(mean(pos(1,n(:,i))), mean(pos(2,n(:,i))), num2str(i), "HorizontalAlignment", "center", "VerticalAlignment","middle", "BackgroundColor", [1,1,1,0.8]);
                end
            end
        end

        % plot geometry
        function plotGeo(obj)
            pdegplot(obj.model, VertexLabels="on", EdgeLabels='on', FaceLabels='on');
        end

        % compute (a_l, b_l, c_l) for the l=1,2,3 basis functions defined
        % on element e
        % the result is index (l, coef), where a is coef=1, b is coef=2, c
        % is coef=3
        function abc = elemABC(obj, e)
            % get node positions
            elems = obj.mesh.Elements;
            n = obj.mesh.Nodes;
            x = n(1,elems(:,e));
            y = n(2,elems(:,e));

            abc = [(x(2)*y(3)-x(3)*y(2)) (y(2)-y(3)) (x(3)-x(2));
                   (x(3)*y(1)-x(1)*y(3)) (y(3)-y(1)) (x(1)-x(3));
                   (x(1)*y(2)-x(2)*y(1)) (y(1)-y(2)) (x(2) - x(1))];
        end

        % compute the area of element e
        function delta_e = areaElem(obj, e)
            abc = obj.elemABC(e);
            delta_e = 0.5*(abc(1,2)*abc(2,3) - abc(2,2)*abc(1,3));
        end

        % evaluate the l-th (l=1,2,or 3) basis function over element e at
        % position (x,y) (which must be inside element e)
        % notice that l is indexed per element, not by global nodal index
        function phi_l = basisElemL(obj, e, l, x, y)
            abc = obj.elemABC(e);
            delta_e = obj.areaElem(e);

            phi_l = (abc(l,1) + abc(l,2)*x + abc(l,3)*y) / 2 / delta_e;
        end

        % get a list of nodes which are subject to dirichlet coditions and
        % their values
        function [dNodes, dVals] = dirichletNodes(obj)
            % outside has potential at 1
            boundaryEdges1 = findNodes(obj.mesh, "region", "Edge", 1:4);
            dNodes1 = unique(boundaryEdges1);
            dVals1 = ones(size(dNodes1));
            % inside has potential at 0
            boundaryEdges2 = findNodes(obj.mesh, "region", "Edge", 5:8);
            dNodes2 = unique(boundaryEdges2);
            dVals2 = zeros(size(dNodes1));
            
            dNodes = [dNodes1 dNodes2].';
            dVals = [dVals1 dVals2].';
        end

        % evaluate K_alpha contributions
        function obj = assembleKalpha(obj)
            elems = obj.mesh.Elements;
            for e=1:length(elems)
                % get nodes of element
                node_indices = elems(:,e);
                delta_e = obj.areaElem(e);
                abc = obj.elemABC(e);
                for l = 1:length(node_indices)
                    for k = 1:length(node_indices)
                        % get corresponding indices in system matrix
                        i = node_indices(l); j = node_indices(k);
                        
                        % inner product of grad of l-th and k-th elements
                        % on e
                        obj.K(i,j) = obj.K(i,j) + 1/(4*delta_e) * (abc(l,2)*abc(k,2) + abc(l,3)*abc(k,3));
                    end
                end
            end
        end

        % correct system matrix and right hand side for dirichlet nodes
        function obj = correctDirichlet(obj)
            [dNodes, dVals] = obj.dirichletNodes();
            for ni=1:length(dNodes)
                n = dNodes(ni);
                % remove the pde enforcement at the boundary and just
                % directly enforce the dirichlet condition instead
                obj.b(n) = dVals(ni);
                obj.K(n,:) = 0;
                obj.K(n,n) = 1;
            end
        end

        % assemble the system matrix
        function obj = assemble(obj)
            % 
            % init empty system matrix
            obj.K = zeros(length(obj.mesh.Nodes));
            obj.b = zeros(length(obj.mesh.Nodes),1);

            % assemble terms in matrix
            obj = obj.assembleKalpha();

            obj = obj.correctDirichlet();

            phi = obj.K \ obj.b;
            pdeplot(obj.mesh, XYData=phi, ColorMap="jet");
        end
    end
end