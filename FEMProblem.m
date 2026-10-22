% 
classdef FEMProblem
    properties
        alpha % pde parameter (alpha)
        beta  % pde parameter
        f % forcing

        geometry % geometry of the problem
        dirEdges % dirichlet BC edeges
        dirVals  % dirichlet BC values
        neuEdges % Neumann BC edges
        neuKappa = {} % Neumann BC values
        neuGamma % Neumann BC values

        model % pde model object
        mesh % pde mesh object
        K % system matrix
        b % right hand side
        phi % solution matrix
    end

    methods

        function obj = FEMProblem(alpha, beta, f)
            obj.alpha = alpha;
            obj.beta = beta;
            obj.f = f;
        end

        function obj = makeRing(obj, geo_a, geo_b)
        % MAKERING  Build a ring and attach it to the model.
        % obj = obj.makeRing(geo_a, geo_b)
        %
        % INPUTS
        % geo_a      inner radius [m]
        % geo_b      outer radius [m]
        %
        % OUTPUTS
        % obj        FEMProblem object with updated obj.model and obj.geometry

            circle_code = 1;
            circle_footer = [];

            % create model
            obj.model = createpde();

            % outer circle (boundary)
            % [code, x-center, y-center, radius, footer]
            B = [circle_code, 0, 0, geo_b, circle_footer]';

            % inner circle (scatterer)
            % [code, x-center, y-center, radius, footer]
            A = [circle_code, 0, 0, geo_a, circle_footer]';

            % TODO: add parameter to mesh and tag (rather than exclude) scatterer (for dielectric case)
            set_names = ['B'; 'A']';
            set_formula = 'B - A';

            % build geometry
            gd = [B A];
            g = decsg(gd, set_formula, set_names);
            obj.geometry = geometryFromEdges(obj.model, g);
        end

        function obj = makeSlice(obj, geo_a, geo_b, swp_angle)
        % MAKESLICE  Build a wedge (slice) and attach it to the model.
        % obj = obj.makeSlice(geo_a, geo_b, swp_angle)
        %
        % INPUTS
        % geo_a      inner radius [m]
        % geo_b      outer radius [m]
        % swp_angle  sweep angle of the slice [rad], 0 < swp_angle <= 2*pi
        %
        % OUTPUTS
        % obj        FEMProblem object with updated obj.model and obj.geometry

            circle_code = 1;
            polygon_code = 2;
            
            % create model
            obj.model = createpde();

            % mask to turn the ring into a wedge
            vertices = 6;
            vert_x = [0; 2*geo_b*cos(linspace(0,swp_angle,5))'];
            vert_y = [0; 2*geo_b*sin(linspace(0,swp_angle,5))'];
            M = [polygon_code; vertices; vert_x; vert_y];
            circle_footer = zeros(1, 2*vertices-2); % adjust the circle footer.

            % outer circle (boundary)
            % [circle, x-center, y-center, radius, footer]
            B = [circle_code, 0, 0, geo_b, circle_footer]';

            % inner circle (scatterer)
            % [circle, x-center, y-center, radius, footer]
            A = [circle_code, 0, 0, geo_a, circle_footer]';

            % TODO: add parameter to mesh and tag (rather than exclude) scatterer (for dielectric case)
            set_names = ['B'; 'A'; 'M']';
            set_formula = '(B - A) * M';

            % build geometry
            gd = [B A M];
            g = decsg(gd, set_formula, set_names);
            obj.geometry = geometryFromEdges(obj.model, g);
        end
        
        function obj = makeRectangles(obj, h, w)
        % MAKERECT  Build two appended rectangles and attach it to the model.
        % obj = obj.makeRect(h, w)
        %
        % INPUTS
        % h          height in y-axis direction [m]
        % w          width in x-axis direction [m]
        %
        % OUTPUTS
        % obj        FEMProblem object with updated obj.model and obj.geometry

            polygon_code = 2;

            % create model
            obj.model = createpde();

            % rectangle (boundary), centered at the origin
            % [polygon, number of vertices, x-coordinates, y-coordinates]
            vertices = 4;
            R1 = [polygon_code; vertices; -w/2; 0;   0;   -w/2; -h/2; -h/2; h/2; h/2];  % left half
            R2 = [polygon_code; vertices;  0;   w/2; w/2;  0;   -h/2; -h/2; h/2; h/2];  % right half

            % build geometry
            g = decsg([R1 R2]);
            obj.geometry = geometryFromEdges(obj.model, g);
        end

        function obj = createMesh(obj, h)
        % CREATEMESH Mesh obj.geometry and attach the mesh to the model.
        % Requires obj.geometry to be non-empty.
        % obj = obj.creatMesh(h)
        %
        % INPUTS
        % h          target element size [m] sets Hmin = h and Hmax = 2*h
        %
        % OUTPUTS
        % obj        FEMproblem object with updated obj.mesh

            if isempty(obj.geometry)
                disp("Meshing Error: Problem geometry needs to be defined before meshing.");
            else
                % mesh after creating the geometry
                obj.mesh = generateMesh(obj.model, Hmax=h*2, Hmin=h, Hgrad=2, GeometricOrder="linear");
            end
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
            
            % figure(); % this adds a figure window we do not need for the
            % manner in which it is called
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

        function obj = setDirichletBC(obj, edges, val)
        % SETDIRICHLETBC Set Dirichlet boundary conditions
        % obj = obj.setDirichletBC(edges, val)
        %
        % INPUTS
        % edges      the edges to which the boundary condition applies
        % val        the value of the boundary condition
        %
        % OUTPUTS
        % obj        FEMproblem object with updated obj.dirEdges and
        %            obj.dirVals

            n = numel(edges);
            obj.dirEdges(end+1:end+n) = edges;
            obj.dirVals(end+1:end+n) = val;
        end

        function obj = setNeumannBC(obj, edges, gammaVal, kappaVal)
        % SETNEUMANNBC Set Neumann boundary conditions
        % obj = obj.setNeumannBC(edges)
        %
        % INPUTS
        % edges      the edges to which the boundary condition applies
        %
        % OUTPUTS
        % obj        FEMproblem object with updated obj.neuEdges and
        %            obj.neuKappa
        
            % if no gamma value specified, default to 0
            if nargin < 3
                gammaVal = 0; 
            end
            
            % if no kappa_n value specified, default to 0
            if nargin < 4
                kappaVal = 0; 
            end

            n = numel(edges);
            obj.neuEdges(end+1:end+n) = edges;
            obj.neuGamma(end+1:end+n) = gammaVal;
            obj.neuKappa(end+1:end+n) = {kappaVal};
        end

        function [dNodes, dVals] = dirichletNodes(obj)
        % DIRICHLETNODES Catalogue the nodes and values corresponding to 
        % edges with dirichlet boundary conditions.
        % obj = obj.dirichletNodes()
        %
        % INPUTS
        % 
        % OUTPUTS
        % dNodes     the nodes with a dirichlet boundary condition
        % dVals      the value of the dirichlet boundary condition
        %            corresponding to each node in dNodes

            % initialize
            dNodes = [];
            dVals  = [];
            for d = 1:numel(obj.dirEdges) % sweep over dirichlet assigned edeges
                nodes  = unique(findNodes(obj.mesh, "region", "Edge", obj.dirEdges(d)));
                dNodes = [dNodes; nodes(:)]; % add all unique nodes
                dVals  = [dVals;  obj.dirVals(d)*ones(numel(nodes),1)]; % add the value
            end
        end

        function [nNodes, nKappa, nGamma] = neumannNodes(obj)
        % NEUMANNNODES Catalog the nodes, and values corresponding to 
        % edges with Neumann boundary conditions.
        % obj = obj.neumannNodes()
        %
        % INPUTS
        % 
        % OUTPUTS
        % nNodes     the nodes with a Neumann boundary condition
        % nKappa     the value of Kappa_n corresponding to each node in nNodes
        % nGamma     the value of Gamma corresponding to each node in nNodes

            % initialize
            nNodes = [];
            nKappa  = [];
            nGamma = [];

            for n = 1:numel(obj.neuEdges) % sweep over neumann assigned edeges
                nodes  = unique(findNodes(obj.mesh, "region", "Edge", obj.neuEdges(n)));
                nNodes = [nNodes; nodes(:)]; % add add unique nodes
                nKappa  = [nKappa;  obj.kappaNode(n, nodes(:))]; % add the evaluated values
                nGamma = [nGamma;  obj.neuGamma(n)*ones(numel(nodes),1)]; % add the value
            end
        end

        function alpha_e = alphaElem(obj, e)
        % ALPHA_E assign alpha values in the scenario of variable or 
        % constant obj.alpha.
        % alpha_e = obj.alphaElem(e)
        %
        % INPUTS
        % e          the element for which obj.alpha is to be
        %            evaluated
        %
        % OUTPUTS
        % alpha_e    if obj.alpha is a function handle this is the function 
        %            evaluated at the centroid of the element cell. If it is
        %            a constant it is set equal to that value.

            if isa(obj.alpha, 'function_handle') 

                elems = obj.mesh.Elements; % elements
                n = obj.mesh.Nodes; % nodes

                alpha_e = zeros(1,size(e,2)); % initialize
                for k = 1:numel(e)
                    node_indices = elems(:,e(k)); % indices for the element

                    % fine the centroid
                    centroidX = (n(1, node_indices(1)) + n(1, node_indices(2)) ...
                        + n(1, node_indices(3))) / 3;
                    centroidY = (n(2, node_indices(1)) + n(2, node_indices(2)) ...
                        + n(2, node_indices(3))) / 3;

                    % evaluate alpha
                    if nargin(obj.alpha) == 2 % cartesian alpha
                        alpha_e(k) = obj.alpha(centroidX, centroidY);
                    else
                        % radial alpha
                        centroid_rho = sqrt(centroidX^2 + centroidY^2);
                        alpha_e(k) = obj.alpha(centroid_rho);
                    end
                end
            else
                alpha_e = obj.alpha * ones(1,size(e,2)); % constant
            end
        end

        function kappa_n = kappaNode(obj, q, nodes)
        % KAPPA_N assign kappa_n values in the scenario of variable or
        % constant obj.neuKappa.
        % kappa_n = obj.kappaNode(q, nodes)
        %
        % INPUTS
        % q          the index of the Neumann edge for which obj.neuKappa 
        %            is to be evaluated
        % nodes      the nodes on that edge for which obj.neuKappa is to be
        %            evaluated
        %
        % OUTPUTS
        % kappa_n    if obj.neuKappa{q} is a function handle this is the function
        %            evaluated at each node. If it is a constant it is set equal
        %            to that value.
            if isa(obj.neuKappa{q}, 'function_handle')

                n = obj.mesh.Nodes; % nodes
                kappa_n = zeros(numel(nodes),1); % initialize

                for k = 1:numel(nodes)

                    % find the node position
                    nodeX = n(1, nodes(k));
                    nodeY = n(2, nodes(k));

                    % evaluate kappa_n
                    if nargin(obj.neuKappa{q}) == 2 % cartesian kappa_n
                        kappa_n(k) = obj.neuKappa{q}(nodeX, nodeY);
                    else
                        % radial kappa_n
                        node_rho = sqrt(nodeX^2 + nodeY^2);
                        kappa_n(k) = obj.neuKappa{q}(node_rho);
                    end
                end
            else
                kappa_n = obj.neuKappa{q} * ones(numel(nodes),1); % constant
            end
        end

        function obj = plotEps(obj)
            [p,e,t] = meshToPet(obj.mesh);
            pdeplot(p, e, t, XYData=obj.alphaElem(1:size(t,2)), XYStyle="flat", ColorMap="jet")
        end

        function obj = assembleKalpha(obj)
        % ASSEMBLEKALPHA Add K_alpha to the K matrix.
        % obj = obj.assembleKalpha()
        %
        % INPUTS
        %
        % OUTPUTS
        % obj        FEMproblem pobject with obj.K additively updated with
        %            the K_alpha term

            elems = obj.mesh.Elements;
            for e=1:size(elems,2)
                % get nodes of element
                node_indices = elems(:,e);
                delta_e = obj.areaElem(e);
                abc = obj.elemABC(e);
                alpha_e = obj.alphaElem(e);

                for l = 1:length(node_indices)
                    for k = 1:length(node_indices)
                        % get corresponding indices in system matrix
                        i = node_indices(l); j = node_indices(k);
                        
                        % inner product of grad of l-th and k-th elements on e

                        % UPDATE: alpha_e is used in place of
                        % obj.alpha to handle cases in which
                        % alpha is varied radially (e.g. test case 3).

                        obj.K(i,j) = obj.K(i,j) ...
                            + alpha_e * (1/(4*delta_e) * (abc(l,2)*abc(k,2) + abc(l,3)*abc(k,3)));
                    end
                end
            end
        end

        function obj = assembleKbeta(obj)
        % ASSEMBLEKBETA Add K_beta to the K matrix.
        % obj = obj.assembleKbeta()
        %
        % INPUTS
        %
        % OUTPUTS
        % obj        FEMproblem pobject with obj.K additively updated with
        %            the K_beta term

            elems = obj.mesh.Elements;
            for e=1:size(elems,2) % for each element

                idx = elems(:,e); % get the nodes from element e
                delta_e = obj.areaElem(e); % area of element e
                        
                % product of l-th and k-th elements on e
                obj.K(idx,idx) = obj.K(idx,idx) ...
                    + obj.beta * delta_e/12 * [2, 1, 1; 1, 2, 1; 1, 1, 2];
            end
        end

        function obj = assembleKgammaBn(obj)
        % ASSEMBLEKGammaBn Add K_gamma to the K matrix and beta_n to the
        % beta vector.
        % obj = obj.assembleKbetaBn()
        %
        % INPUTS
        %
        % OUTPUTS
        % obj        FEMproblem with obj.K additively updated with the 
        %            K_gamma term and obj.b updated with the beta_n term.

            elems = obj.mesh.Elements;
            [nNodes, nKappa, nGamma] = obj.neumannNodes();

            for e=1:size(elems,2) % for each element
                node_indices = elems(:,e);
                
                for l = 1:length(node_indices) % for each node
                    for k = 1:length(node_indices) % for each node
                        i = node_indices(l); j = node_indices(k);

                        if i == j
                            continue
                        end

                        iN = find(nNodes == i, 1, 'first');
                        jN = find(nNodes == j, 1, 'first');

                        if (numel(iN) == 0) || (numel(jN) == 0)
                            continue
                        end

                        gammaN = nGamma(iN);
                        kappaN = nKappa(iN);

                        ijLen = sqrt((obj.mesh.Nodes(1,i) - obj.mesh.Nodes(1,j))^2 + (obj.mesh.Nodes(2,i) - obj.mesh.Nodes(2,j))^2);

                        % add K contribution
                        obj.K(i,j) = obj.K(i,j) + ijLen * 1/6 * gammaN;
                        obj.K(i,i) = obj.K(i,i) + ijLen * 1/3 * gammaN;

                        % add contributions to beta_n
                        obj.b(i) = obj.b(i) + 0.5 * kappaN * ijLen;
                    end
                end
            end
        end

        function obj = assemblebF(obj)
        % ASSEMBLEBF Add b_F to the b matrix.
        % obj = obj.assemblebF()
        %
        % INPUTS
        % 
        % OUTPUTS
        % obj        FEMproblem pobject with obj.b additively updated with
        %            the b_F term

            elems = obj.mesh.Elements;
            for e=1:size(elems,2) % for each element

                idx = elems(:,e); % get the nodes from element e
                delta_e = obj.areaElem(e); % area of element e
                        
                % product of l-th and k-th elements on e
                obj.b(idx) = obj.b(idx) + obj.f * delta_e/3;
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
            obj = obj.assembleKbeta();
            obj = obj.assembleKgammaBn();
            obj = obj.assemblebF();
            obj = obj.correctDirichlet();

            obj.phi = obj.K \ obj.b;
            pdeplot(obj.mesh, XYData=abs(obj.phi), ColorMap="jet");
        end

        function Cap = calculateCap(obj)
            % Calculate capacitance:
            elems = obj.mesh.Elements;
            [gradX, gradY] = pdegrad(obj.mesh.Nodes, elems, obj.phi); % calculate gradients

            W = 0; % electric energy
            for e = 1:size(elems,2)
                alpha_e = obj.alphaElem(e); % update obj.alpha for element e
                W = W + alpha_e * obj.areaElem(e) * (gradX(e)^2 + gradY(e)^2);
            end

            % Capacitance (calculate and save)
            Cap = W * 8.85 * 1e-12; % F/m per Volt  
        end
    end
end