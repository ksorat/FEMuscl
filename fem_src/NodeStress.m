function sig_n = NodeStress(N,eltype,numels,numnodes,sig_v,L,J,LM,w)

global dim Elements Nodes Parts;
mM = zeros(dim*length(Elements));
P = zeros(dim*length(Elements),1);


for i = 1:length(Elements)

    pe = zeros((dim-1)*3*length(Elements{i}),1);
    m = zeros((dim-1)*3*length(Elements{i}));
    Nm = Parts{Elements{i}.part}.elForm.StressMatrix;
    Ns = Parts{Elements{i}.part}.elForm.Ns;
    
    for j = 1:length(Nm)
        [~,~,J] = Elements{i}.BMatrix(j);
        pe = pe + Ns{j}'*Elements{i}.VoigtStress(j); % Calculate element load vector
        m = m + Nm{j}*J; % Calculate element "mass" matrix
    end
    
    for i = 1:length(pe)
               
       P = P + pe(i);

    end % scatter nodal forces into global vector
    
    P = P + L{i}'*pe; % Assemble Global Load Vector
    mM = mM + L{i}'*m*L{i}; 
    
end

sig_hat = mM\P; % Solve for nodal stresses using L2
sig_n = stresscomps(numnodes,sig_hat); % separate out stresses for each node